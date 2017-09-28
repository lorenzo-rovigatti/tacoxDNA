# tacoxDNA

tacoxDNA (Tools and Converters for oxDNA) is a collection of tools developed to help [oxDNA](http://dna.physics.ox.ac.uk/) users. The sections that follow introduce the tools and their usage.

## Generator for twisted and knotted configurations

## oxDNA-to-LAMMPS converter

The '''oxDNA-to-LAMMPS.py''' script takes two mandatory arguments and outputs a single file.

#### Arguments
* An oxDNA topology file
* An oxDNA configuration file

#### Output
* A file containing the list of nucleotide positions, quaternions, velocities, angular velocities and bonds which can be used as a restart file in LAMMPS. The name of the file is just the oxDNA configuration filename, prefixed with "LAMMPS_"

## LAMMPS-to-oxDNA converter

The '''LAMMPS-to-oxDNA.py''' script takes two mandatory arguments and outputs two files.

#### Arguments
* A LAMMPS input file, containing the topology of the configuration in the form of a list of bonds
* A LAMMPS output file having [...]

#### Output
* An oxDNA topology file (named by suffixing the LAMMPS output file with ".top")
* An oxDNA configuration file (named by suffixing the LAMMPS output file with ".conf")

## oxDNA-to-all-atom converter

The '''oxDNA-to-all-atom.py" script takes two mandatory arguments and outputs a single file.

#### Arguments
* An oxDNA topology file
* An oxDNA configuration file

#### Output
* A [PDB](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html) file containing the positions of all atoms composing the strand(s) contained in the oxDNA configuration file 

## Twist and writhe analyzer

## Acknowledgements

* Some of the files in the src/libs folder have been adapted from oxDNA source code
