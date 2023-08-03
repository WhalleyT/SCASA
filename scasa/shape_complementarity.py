import math
import random

import numpy as np
import scipy.spatial

from scipy.spatial import cKDTree
from dataclasses import dataclass
from itertools import compress
from scipy.spatial import Delaunay
from sklearn.decomposition import PCA


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
        est = scipy.spatial.ConvexHull(c)
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
        v1 = p2 - p1
        v2 = p3 - p1
        v3 = np.cross(v1, v2)
        return v3 / np.linalg.norm(v3)

    def PCA_unit_vector(array, pca=PCA(n_components=3)):
        pca.fit(array)
        eigenvalues = pca.explained_variance_
        return pca.components_[np.argmin(eigenvalues)]

    def calculate_normal(self, coordinate, mesh):
        tree = cKDTree(mesh)
        # Get indices and distances:
        dist, ind = tree.query(coordinate, k=3)
        combinations = mesh[ind]

        #calculate normals using cross product
        normals = self.calc_cross(combinations[:,0], combinations[:,1], combinations[:,2])
        print(normals)
        normals = self.calc_cross(combinations[0], combinations[1], combinations[2])
        print(normals)
        normals = self.PCA_unit_vector(combinations)
        print(normals)


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

        #for set of random points find nearest point
        sc_a, sc_b = [], []
        # todo calculate normals for neighbouring points

        for coordinate_complex_1 in points_c1:
            coordinate_complex_2 = self.find_nearest_neighbour(coordinate_complex_1, points_c2)

            normal_coord_complex_1 = self.calculate_normal(coordinate_complex_1, points_c1)
            normal_coord_complex_2 = self.calculate_normal(coordinate_complex_2, points_c2)


        #sc = (np.median(np.array(sc_a)) + np.median(np.array(sc_b))) / 2
        #print(sc)



    def get_column(self, param):
        pass

