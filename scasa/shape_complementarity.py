import math
import random

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import plotly.figure_factory as ff
import plotly.express as px

from scipy.spatial import cKDTree, Delaunay, ConvexHull
from sklearn.decomposition import PCA
from dataclasses import dataclass
from itertools import compress


@dataclass
class PDBCoords:
    """
    Dataclass containing PDB coordinates, amino acids, atoms and residue numbers
    """
    coords: np.array
    amino_acids: list
    atoms: list
    residues: list


class ShapeComplementarity:
    def __init__(self, arg):
        self.verbose = None
        self.distance = None
        self.density = None
        self.complex_1 = None
        self.complex_2 = None
        self.arg = arg
        super().__init__(self.arg)

    def convert_1d_array(self, arr):
        """
        return a numpy array from a list
        :param arr: list
        :return: np.array
        """
        return np.array(arr, dtype=float)

    def create_interface(self):
        """
        create an interface of the two objects, i.e. convert them into two separate
        numpy arrays
        :return: two 3xn numpy arrays
        """

        if self.verbose:
            print("Getting interfaces")

        residues = self.get_column("RESIDUE_NUM")
        chain = self.get_column("CHAIN")
        amino_acids = self.get_column("RESIDUE_SEQID")
        atoms = self.get_column("ATOM_NAME")

        x_coord = self.convert_1d_array(self.get_column("X_COORD"))
        y_coord = self.convert_1d_array(self.get_column("Y_COORD"))
        z_coord = self.convert_1d_array(self.get_column("Z_COORD"))

        complex_1_created = False
        complex_2_created = False

        complex_1_residues = []
        complex_2_residues = []

        complex_1_aa = []
        complex_2_aa = []

        complex_1_at = []
        complex_2_at = []
        for x, y, z, c, r, aa, a in zip(x_coord, y_coord, z_coord, chain,
                                        residues, amino_acids, atoms):
            coord = np.array((x, y, z), "f")
            if c in self.complex_1:
                complex_1_residues.append(r)
                complex_1_aa.append(aa)
                complex_1_at.append(a)
                if complex_1_created:
                    complex_1_coords = np.append(complex_1_coords, coord)
                else:
                    complex_1_coords = coord
                    complex_1_created = True
            if c in self.complex_2:
                complex_2_residues.append(r)
                complex_2_aa.append(aa)
                complex_2_at.append(a)
                if complex_2_created:
                    complex_2_coords = np.append(complex_2_coords, coord)
                else:
                    complex_2_coords = coord
                    complex_2_created = True

        complex_1_coords = np.reshape(complex_1_coords, (-1, 3))
        complex_2_coords = np.reshape(complex_2_coords, (-1, 3))

        c1 = PDBCoords(coords=complex_1_coords, amino_acids=complex_1_aa,
                       atoms=complex_1_at, residues=complex_1_residues)
        c2 = PDBCoords(coords=complex_2_coords, amino_acids=complex_2_aa,
                       atoms=complex_2_at, residues=complex_2_residues)
        return c1, c2

    def filter_interface(self, a, b, r):
        """
        Return atoms only within n Angstroms
        :return:
        """

        tree = cKDTree(a.coords)
        mask = np.zeros(len(a.coords), dtype=bool)

        indices = []
        for b in b.coords:
            i = tree.query_ball_point(b, r)
            indices += i

        indices = list(set(indices))
        mask[indices] = True

        # Extract values not among the nearest neighbors
        c = a.coords[mask]
        aa = list(compress(a.amino_acids, list(mask)))
        at = list(compress(a.atoms, list(mask)))
        res = list(compress(a.residues, list(mask)))

        return PDBCoords(coords=c, amino_acids=aa, atoms=at, residues=res)

    def estimate_surface_area(self, c):
        """
        estimate surface area of numpy array using Convex Hull
        :param c: numpy array (in this case 3D
        :return: area, float
        """
        est = ConvexHull(c)
        if self.verbose:
            print(f"Estimated area of complex 1's face is {est.area:.2f}Å\N{SUPERSCRIPT THREE}")
        return est.area

    def create_polygon(self, points):
        """
        Converts 3D array of a surface into a Delaunay triangulation. returns simplexes
        :param points: np.array 
        :return: scipy.spatial.Delaunay.simplices.
        array of indexes pointing to coordinates that constitute a triangle
        """
        x = points[:, 0]
        y = points[:, 1]
        z = points[:, 2]
        # todo check if x,y,z needs to be specificed so that interface is in correct orientation
        two_d = np.vstack([x, y]).T
        tri = Delaunay(two_d)
        simplices = tri.simplices
        return simplices


    def point_inside_triangle(self, v1, v2, v3):
        """
        randomly generate point inside of triangle where:
        P = (1−a−−√)v1 + (a−−√(1−b))v2 + (ba−−√)
        :param v1: vertex 1 of triangle, a numpy array of 3 coordinates
        :param v2: vertex 2
        :param v3: vertex 3
        :return: np.array of 3 coordinates. float.
        """
        #P = (1−a−−√)v1 + (a−−√(1−b))v2 + (ba−−√)

        a = random.random()
        b = random.random()
        p = (1- math.sqrt(a)) * v1 + (math.sqrt(a) * (1 - b)) * v2 \
              + (b * math.sqrt(a)) * v3
        return p

    def random_points(self, coords, simplexes, n_samp):
        """
        Generate random points for n_samples using a triangluation of a mesh
        :param coords: np.array 3xn
        :param simplexes: np.array 3xn
        :param n_samp: number of random samples
        :return: 3xn_samp np.array
        """
        # get random indexes of size n_samp
        indices = np.random.choice(coords.shape[0], n_samp)
        #create a simplex matrix from the sample
        simp = simplexes[indices]

        #for a given simplex, sample a random point
        random_points = []
        for s in simp:
            #extract 3 points of a triangle
            p1 = coords[s[0]]
            p2 = coords[s[1]]
            p3 = coords[s[2]]
            p = self.point_inside_triangle(p1, p2, p3)
            random_points.append(p)

        return np.array(random_points)

    def surface_complementarity_function(self, x_a, xd_a, w):
        """
        compute S(A->B)(x A) = (n A.n’ A) exp [-w(|x A - x’ A |) 2]
        :param x_a: point on the interface of surface A
        :param xd_a: nearest point to x_a on surface B
        :param w: weighting factor
        :return: float
        """
        return np.dot(x_a, xd_a) * np.exp(-w*(np.linalg.norm(x_a -xd_a)))

    def find_nearest_neighbour(self, coord, set_of_coords):
        """

        :param coord: single x,y,z coord as np.array
        :param set_of_coords:  array of x,y,z coord
        :return: coord, np.array of nearest coordinate
        """
        tree = cKDTree(set_of_coords)
        nearest_coord = tree.query(coord)

        return set_of_coords[nearest_coord[1]]

    def calc_cross(self, p1, p2, p3):
        """
        calculate cross product of three points
        :param p1: point 1 np.array
        :param p2: point 2 np.array
        :param p3: point 2 np.array
        :return: cross product np.array
        """
        v1 = p2 - p1
        v2 = p3 - p1
        v3 = np.cross(v1, v2)
        return v3 / np.linalg.norm(v3)

    def calculate_normal(self, coordinate, mesh):
        """
        Calculate surface normal from a point using its two neighbours.
        I have used PCA as I can arbitarly adjust n neighbours.
        :param coordinate: coordinate of interest
        :param mesh: mesh that the coordinate of interest is from.
        :return: np.array normals, [x,y,z]
        """
        #todo check if the positive/negativity of this matters
        #create
        tree = cKDTree(mesh)
        # Get indices and distances. since dist=0 we will return our point of interest so no need to add
        dist, ind = tree.query(coordinate, k=5)
        #get np array of sub indexes
        combinations = mesh[ind]

        #calculate cross product
        pca = PCA(n_components=3)
        pca.fit(combinations)
        eigenvalues = pca.explained_variance_
        return pca.components_[np.argmin(eigenvalues)]

    def plot_sc(self, sc_complex_1, sc_complex_2):
        complex_1 = pd.DataFrame(sc_complex_1, columns=["SC_function"])
        complex_1["Complex"] = "Complex 1"
        complex_2 = pd.DataFrame(sc_complex_2, columns=["SC_function"])
        complex_2["Complex"] = "Complex 2"

        complexes = pd.concat([complex_1, complex_2])

        sns.histplot(data=complexes, x="SC_function", hue="Complex")
        plt.show()


    def calculate_sc(self, points_c1, points_c2, weight):
        """
        Calculates the surface complementarity function for every point in an array of points for a single complex.
        That is SC(points_c1 -> points_c2)
        :param points_c1: complex one, the part of the protein you are calculating for. np.array (3xn)
        :param points_c2: complex two, the part of the protein you are calculating against. np.array (3xn)
        :param weight: weighting parameter
        :return: list of SC scores (float)
        """
        sc_array = []

        for coordinate_complex_1 in points_c1:
            coordinate_complex_2 = self.find_nearest_neighbour(coordinate_complex_1, points_c2)

            normal_coord_complex_1 = self.calculate_normal(coordinate_complex_1, points_c1)
            normal_coord_complex_2 = self.calculate_normal(coordinate_complex_2, points_c2)

            distance_function = np.exp(-(np.linalg.norm(coordinate_complex_1-coordinate_complex_2)**2) * weight)

            dot_prod = np.dot(normal_coord_complex_1, normal_coord_complex_2)
            sc = dot_prod * distance_function
            print(sc, dot_prod, distance_function)
            sc_array.append(sc)
        return sc_array


    def plot_combined_mesh(self, mesh1, mesh2, coords1, coords2):
        x_1 = coords1[:, 0]
        y_1 = coords1[:, 1]
        z_1 = coords1[:, 2]

        x_2 = coords2[:, 0]
        y_2 = coords2[:, 1]
        z_2 = coords2[:, 2]

        x = np.concatenate([x_1, x_2])
        y = np.concatenate([y_1, y_2])
        z = np.concatenate([z_1, z_2])

        mesh = np.concatenate([mesh1, mesh2])

        colours = ["#EF553B"] * len(x_1) + ["#00CC96"] * len(x_2)

        fig = ff.create_trisurf(x=x,y=y,z=z, simplices=mesh, colormap=colours, title="Surface of  Complex 1 and 2")
        return fig

    def plot_single_mesh(self, mesh, coords, title):
        x = coords[:, 0]
        y = coords[:, 1]
        z = coords[:, 2]

        fig = ff.create_trisurf(x=x,y=y,z=z, simplices=mesh, title=title)
        return fig

    def plot_atoms(self, c1, c2, title):
        x_1 = c1[:, 0]
        y_1 = c1[:, 1]
        z_1 = c1[:, 2]

        x_2 = c2[:, 0]
        y_2 = c2[:, 1]
        z_2 = c2[:, 2]

        x = np.concatenate([x_1, x_2])
        y = np.concatenate([y_1, y_2])
        z = np.concatenate([z_1, z_2])

        colours = ["Complex 1"] * len(x_1) + ["Complex 2"] * len(x_2)

        df = pd.DataFrame(list(zip(x,y,z, colours)), columns=["X", "Y", "Z", "Complex"])

        fig = px.scatter_3d(df, x='X', y='Y', z='Z',
                            color='Complex', title=title)
        return  fig



    def sc(self):
        complex1, complex2 = self.create_interface()
        complex1 = self.filter_interface(complex1, complex2, self.distance)
        complex2 = self.filter_interface(complex2, complex1, self.distance)

        if self.verbose:
            print(f"Complex 1 contains {len(complex1.residues)} atoms "
                  f"within {self.distance} Angstroms of Complex 2")
            print("Complex 2 contains %i atoms within %i Angstroms of Complex 1" 
                  % (len(complex2.residues), self.distance))

        area_1 = self.estimate_surface_area(complex1.coords)
        area_2 = self.estimate_surface_area(complex2.coords)

        simplices_c1 = self.create_polygon(complex1.coords)
        simplices_c2 = self.create_polygon(complex2.coords)

        n_dots_1 = round(self.density * area_1)
        n_dots_2 = round(self.density * area_2)

        points_c1 = self.random_points(complex1.coords, simplices_c1, n_dots_1)
        points_c2 = self.random_points(complex2.coords, simplices_c2, n_dots_2)

        if self.verbose:
            print("Calculating SC for both complexes")

        sc_complex_1 = self.calculate_sc(points_c1, points_c2, self.weight)
        sc_complex_2 = self.calculate_sc(points_c2, points_c1, self.weight)

        sc_score = (np.median(sc_complex_1) + np.median(sc_complex_2) ) / 2
        if self.verbose:
            print(f"SC = {sc_score:.2f}")

        if self.plot:
            p1 = self.plot_atoms(complex1.coords, complex2.coords, "Complex 1 and 2 selected atoms")
            p2 = self.plot_combined_mesh(simplices_c1, simplices_c2, complex1.coords, complex2.coords)
            p3 = self.plot_single_mesh(simplices_c1, complex1.coords, "Surface of  Complex 1")
            p4 = self.plot_single_mesh(simplices_c2, complex2.coords, "Surface of  Complex 2")
            p5 = self.plot_atoms(points_c1, points_c2, "Complex 1 and 2 randomly sampled atoms")

            for p in [p1, p2, p3, p4, p5]:
                p.show()

    def get_column(self, param):
        pass

