# tacox<img src="logo.png" width="35">NA

tacoxDNA (Tools and Converters for oxDNA) is a collection of tools initially developed to help [oxDNA](http://dna.physics.ox.ac.uk/) users. However, it will soon be expanded so as to support additional models. If you use tacoxDNA, please consider citing the following article:

A. Suma, E. Poppleton, M. Matthies, P. Šulc, F. Romano, A. A. Louis, J. P. K. Doye, C. Micheletti and L. Rovigatti, ["TacoxDNA: A user‐friendly web server for simulations of complex DNA structures, from single strands to origami"](https://doi.org/10.1002/jcc.26029), *J. Comput. Chem.* **40**, 2586 (2019)

The sections that follow introduce the tools and their usage.

* [Generator for twisted and knotted configurations](#generator-for-twisted-and-knotted-configurations)
* [oxDNA-to-LAMMPS converter](#oxdna-to-lammps-converter)
* [LAMMPS-to-oxDNA converter](#lammps-to-oxdna-converter)
* [oxDNA-to-PDB converter](#oxdna-to-pdb-converter)
* [PDB-to-oxDNA converter](#pdb-to-oxdna-converter)
* [cadnano-to-oxDNA converter](#cadnano-to-oxdna-converter)
* [CanDo-to-oxDNA converter](#cando-to-oxdna-converter)
* [Tiamat-to-oxDNA converter](#tiamat-to-oxdna-converter)
* [vHelix-to-oxDNA converter](#vhelix-to-oxdna-converter)
* [rpoly-to-oxDNA converter](#rpoly-to-oxdna-converter)

---

## Generator for twisted and knotted configurations

The `XYZ_oxDNA.py` script generates an oxDNA topology/configuration pair from a file containing a list of coordinates that defines a centerline.

### Mandatory arguments
* A centerline file containing a list of coordinates in the format x y z

### Optional arguments
* `-c\--closed`
the last bead is connected to the first bead (default)
* `-o\--open`
the last bead is not connected to the first bead
* `-h\--help`
print usage
* `-d\--dsDNA`
the chain is clad with a double-stranded DNA (default)
* `-s\--ssDNA`
the chain is clad with a single-stranded DNA
* `-n\--nicked`
optional argument when -d option is used. One of the two strands of double strand DNA is nicked (not circularized) (not set by default)
* `-p\--supercoiling=SUPERCOILING_DENSITY`
supercoiling density percentage (defaults to 0, with an equilibrium pitch of 10.5 imposed)
* `-w\--writhe=WRITHE_AMOUNT`
topological target writhe to superimpose. When the supercoiling density is zero, this number corresponds to the average writhe. Useful for knots which have an average writhe different from 0 (it defaults to 0)
* `-e\--seed=RNG_SEED`
random seed for DNA sequence (defaults to a random value)
* `-q\--sequence=SEQUENCE`
text file containing a valid DNA sequence (*e.g.* ATCTGA). The length of the sequence should correspond to the number of points in the coordinate file. If not specified, the sequence will be chosen randomly

### Output
* An oxDNA topology file (named by suffixing the centerline file with ".top")
* An oxDNA configuration file (named by suffixing the centerline file with ".oxdna")

---

## oxDNA-to-LAMMPS converter

The `oxDNA_LAMMPS.py` script takes two mandatory arguments and outputs a single file that can be used as input to run oxDNA simulations [USER-CGDNA](https://lammps.sandia.gov/doc/Packages_details.html#pkg-user-cgdna) LAMMPS package.

### Arguments
* An oxDNA topology file
* An oxDNA configuration file

### Output
* A file containing the list of nucleotide positions, quaternions, velocities, angular velocities and bonds which can be used as a start file in LAMMPS. The name of the file is just the oxDNA configuration filename, prefixed with "LAMMPS_"

---

## LAMMPS-to-oxDNA converter

The `LAMMPS_oxDNA.py` script takes one mandatory argument, one optional argument and outputs two files.

### Arguments
* A LAMMPS data file, containing the nucleotide positions, quaternions, velocities, angular momenta and the bond list
* Optionally, a LAMMPS trajectory file with atom data for each time step after the data file name

### Output
* An oxDNA topology file (named by suffixing the LAMMPS output file with ".top")
* An oxDNA configuration file (named by suffixing the LAMMPS output file with ".oxdna"). If a trajectory file is also processed, the configuration file will contain the full trajectory data in native oxDNA format.

---

## oxDNA-to-PDB converter

The `oxDNA_PDB.py` script takes three mandatory arguments and outputs a single file. Since oxDNA bases has no one-to-one explicit mapping to all-atom representations, the converted structure will most likely require some sort of relaxation procedure before being used as input for all-atom simulation packages. Moreover, when using this script the following points should be taken into account:

* the phosphate groups of nucleobases at the 5' end of each strand are removed
* a hydrogen is added to each of the 3' and 5' end nucleobases (HO3' and HO5', respectively)

### Mandatory arguments
* An oxDNA topology file
* An oxDNA configuration file
* The direction according to which the nucleotides are to be listed in the PDB file. It should be either 35 (for 3' -> 5') or 53 (for 5' -> 3'). Most of the all-atoms tools (*e.g.* GROMACS) assume the 5' -> 3' order.

### Optional arguments
* `-H, --hydrogens=[true|false]` 
if true, include hydrogen atoms in the PDB file (defaults to true)
* `-u\--uniform-residue-names`
drop the `3` and `5` suffixes from the names of residues that are placed at the strands' ends. It increases the compatibility with some tools (*e.g.* Chimera and Molecular Maya)
* `-o\--one-file-per-strand`
print one PDB file for each strand

### Output
* A [PDB](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html) file containing the positions of all atoms composing the strand(s) contained in the oxDNA configuration file

---

## PDB-to-oxDNA converter

The `PDB_oxDNA.py` script takes two mandatory arguments. Given the sometimes messy nature of PDB files, the script makes some choices during the parsing of the input file. In particular, note the following points:

* if the PDB file contains more than one MODEL, only the first one will be converted unless the `-m/--models-as-strands` option is given
* if the PDB file contains alternate locations for some (or all) of the atoms, only those marked with either "1" or "A" will be considered. If the PDB file uses a different notation, the script may fail or crash
* sometimes, sugar atoms are marked with asterisks (\*) instead of single quotes ('). In these cases the converter replaces the former with the latter and moves on

### Mandatory arguments
* The input PDB file
* The direction according to which the nucleotides are listed in the PDB file. It should be either 35 (for 3' -> 5') or 53 (for 5' -> 3').

### Optional arguments
* `-m, --models-as-strands` 
Treat different models as different strands

### Output
* An oxDNA topology file (named by suffixing the PDB file with ".top")
* An oxDNA configuration file (named by suffixing the PDB file with ".oxdna")

---

## cadnano-to-oxDNA converter

The `cadnano_oxDNA.py` script converts [cadnano](https://cadnano.org/) files into oxDNA configurations. Optionally, it can also output [.oxview files](https://sulcgroup.github.io/oxdna-viewer/) which will contain additional information about base pairing, custom colors and clustered domains.

Note that the script **does not** support scaffold-less input files. It takes two mandatory arguments.

### Mandatory arguments
* The input cadnano file (in json format)
* The lattice type the file was designed with. It should be either *sq* (square) or *he* (hexagonal)

### Optional arguments
* `-e\--seed=RNG_SEED`
random seed for DNA sequence (defaults to a random value)
* `-b\--box=VALUE`
the length of the box side (in oxDNA simulation units) where the system will be placed
* `-q\--sequence=SEQUENCE`
text file containing a valid DNA sequence (*e.g.* ATCTGA). If not specified, the sequence will be chosen randomly
* `-p\--print-virt2-nuc`
print the `virt2nuc` file that can be used by the oxDNA's `origami_utils.py` script to convert between cadnano and oxDNA nucleotide indexes
* `-o\--print-oxview`
print a `.oxview` file that can be opened and edited in [oxView](https://sulcgroup.github.io/oxdna-viewer/). Using this option will allow you to keep additional design information not included in the oxDNA files.

### Output
* An oxDNA topology file (named by suffixing the cadnano file with ".top")
* An oxDNA configuration file (named by suffixing the cadnano file with ".oxdna")

---

## CanDo-to-oxDNA converter

The `CanDo_oxDNA.py` script converts [CanDo](https://cando-dna-origami.org/) files into oxDNA configurations. It takes one mandatory argument.

### Mandatory arguments
* The input CanDo file

### Optional arguments
* `-b\--box=VALUE`
the length of the box side (in oxDNA simulation units) where the system will be placed (defaults to 100)
* `-f\--print-force-file`
also print a file containing the specifics for a oxDNA-compatible set of external forces that may be useful to relax the system  

### Output
* An oxDNA topology file (named by suffixing the CanDo file with ".top")
* An oxDNA configuration file (named by suffixing the CanDo file with ".oxdna")

---

## Tiamat-to-oxDNA converter

The `Tiamat_oxDNA.py` script converts [Tiamat](http://yanlab.asu.edu/Resources.html) files into oxDNA configurations. It takes one mandatory argument.

### Mandatory arguments
* The input Tiamat file (in json format)

### Optional arguments
* `-m\--molecule=[DNA|RNA]`
the type of molecule contained in the input file (defaults to DNA)
* `-t\--tiamat-version=[1|2]`
the Tiamat version the input file was generated with. If you are not sure, it's probably 2, the default value
* `-d\--default-base=[A|C|G|T|R|integer]`
some of the bases generated by Tiamat have no associated type. By default, these bases are assigned a random type (either A, C, G or T). By setting this option the user can assign to these bases the same type. Since oxDNA can also use integer numbers as types, these are also supported here
* `-f\--print-force-file`
also print a file containing the specifics for a oxDNA-compatible set of external forces that may be useful to relax the system

### Output
* An oxDNA topology file (named by suffixing the Tiamat file with ".top")
* An oxDNA configuration file (named by suffixing the Tiamat file with ".oxdna")

---

## vHelix-to-oxDNA converter

The `vHelix_oxDNA.py` script converts [vHelix](http://www.vhelix.net/) files into oxDNA configurations. vhelix files are Maya files stored in the MA format.  It takes one mandatory argument.

### Mandatory arguments
* The input vHelix file

### Optional arguments
* `-b\--box=VALUE`
the length of the box side (in oxDNA simulation units) where the system will be placed (defaults to 100)
* `-e\--seed=RNG_SEED`
random seed (defaults to a random value). Random vectors are used whenever the input configuration contains deleted nucleotides.

### Output
* An oxDNA topology file (named by suffixing the vHelix file with ".top")
* An oxDNA configuration file (named by suffixing the vHelix file with ".oxdna")

---

## rpoly-to-oxDNA converter

The `rpoly_oxDNA.py` script converts routed polyhedra (rpoly) files containing wireframe DNA origami structures automatically generated using the BSCOR package(http://www.vhelix.net/) into oxDNA configurations. It takes one mandatory argument.

### Mandatory arguments
* The input rpoly file

### Optional arguments
* `-e\--seed=RNG_SEED`
random seed for DNA sequence (defaults to a random value)

### Output
* An oxDNA topology file (named by suffixing the rpoly file with ".top")
* An oxDNA configuration file (named by suffixing the rpoly file with ".oxdna")

---

## Testing

tacoxDNA contains a very simple testing suite to verify the working status of the scripts. The `tests` directory contains a directory for each script. Within each directory there is a `run.sh` bash script that performs one or more tests on the specific script. Execute `run_all.sh` to run all tests and get a summary of the results. 

## Acknowledgements

* Some of the code has been adapted from the [oxDNA](http://dna.physics.ox.ac.uk/) source
* The vHelix-to-oxDNA converter was provided by Erik Benson
* The source of code of the [pyquaternion lib](https://github.com/KieranWynn/pyquaternion) is included in the source tree
