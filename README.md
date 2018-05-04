# tacoxDNA

tacoxDNA (Tools and Converters for oxDNA) is a collection of tools developed to help [oxDNA](http://dna.physics.ox.ac.uk/) users. The sections that follow introduce the tools and their usage.

## Generator for twisted and knotted configurations

The `centerline-to-oxDNA.py` script generates an oxDNA topology/configuration pair from a file containing a list of coordinates that defines a centerline.

### Mandatory arguments
* A centerline file containing a list of coordinates in the format x y z

### Optional arguments
* `-c\--closed`
the last bead (base) is connected to the first bead (base). Default is True
* `-o\--open`
the last bead (base) is not connected to the first bead (base). Default is False
* `-h\--help`
print usage
* `-d\--dsDNA`
the chain is clad with a double strand DNA. Default is True
* `-s\--ssDNA`
the chain is clad with a single strand DNA. Default is False
* `-n\--nicked`
optional argument when -d option is used. One of the two strands of double strand DNA is nicked (not circularized). Default is False
* `-p\--supercoiling=SUPERCOILING_DENSITY`
supercoiling density percentage. Default is 0, with an equilibrium pitch of 10.5 imposed
* `-w\--writhe=WRITHE_AMOUNT`
additional writhe to superimpose (useful for knots which have an average writhe different from 0). Default value is 0.
* `-e\--seed=RNG_SEED`
random seed.
* `-q\--sequence=SEQUENCE`
text file containing sequence, ex ATCTGA. Default is a random sequence. The length of the sequence should correspond to the number of points in the coordinate file

### Output
* An oxDNA topology file (named by suffixing the centerline file with ".top")
* An oxDNA configuration file (named by suffixing the centerline file with ".conf")

## oxDNA-to-LAMMPS converter

The `oxDNA-to-LAMMPS.py` script takes two mandatory arguments and outputs a single file.

### Arguments
* An oxDNA topology file
* An oxDNA configuration file

### Output
* A file containing the list of nucleotide positions, quaternions, velocities, angular velocities and bonds which can be used as a start file in LAMMPS. The name of the file is just the oxDNA configuration filename, prefixed with "LAMMPS_"

## LAMMPS-to-oxDNA converter

The `LAMMPS-to-oxDNA.py` script takes one mandatory argument and outputs two files.

### Arguments
* A LAMMPS input start file, containing the nucleotide positions, quaternions, velocities, angular velocities and bond list 

### Output
* An oxDNA topology file (named by suffixing the LAMMPS output file with ".top")
* An oxDNA configuration file (named by suffixing the LAMMPS output file with ".conf")

## oxDNA-to-PDB converter

The `oxDNA_PDB.py` script takes three mandatory arguments and outputs a single file.

### Mandatory arguments
* An oxDNA topology file
* An oxDNA configuration file
* The direction according to which the nucleotides are to be listed in the PDB file. It should be either "35" (for 3' -> 5') or 53 (for 5' -> 3').

### Optional arguments
* `-H, --hydrogens=[true|false]` if true, include hydrogen atoms in the PDB file (defaults to true)

### Output
* A [PDB](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html) file containing the positions of all atoms composing the strand(s) contained in the oxDNA configuration file

## PDB-to-oxDNA converter

The `PDB_oxDNA.py` script takes two mandatory arguments.

### Mandatory arguments
* The input PDB file
* The direction according to which the nucleotides are listed in the PDB file. It should be either "35" (for 3' -> 5') or 53 (for 5' -> 3').

### Output
* An oxDNA topology file (named by suffixing the PDB file with ".top")
* An oxDNA configuration file (named by suffixing the PDB file with ".conf")

## Acknowledgements

* Some of the code has been adapted from the [oxDNA](http://dna.physics.ox.ac.uk/) source
