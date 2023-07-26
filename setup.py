import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
     name='scasa',
     version='0.1',
     scripts=['scripts/SCASA'],
     author="Tom Whalley",
     author_email="whalleyt@cardiff.ac.uk",
     description="Python package for the calculation of buried and available surface"
                 " area and shape complementarity of PDB files",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/whalleyt/scasa",
     packages=['scasa'],
    package_dir={"scasa": "scasa"},
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
 )