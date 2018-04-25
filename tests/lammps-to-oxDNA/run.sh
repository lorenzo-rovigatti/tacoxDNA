#!/bin/bash

CORRECT_OUTPUT="correct_output.dat"
CORRECT_TOP="correct_output.top"
OUTPUT="init_lammps.dat.oxdna"
OUTPUT_TOP="init_lammps.dat.top"

if [ ! -s init_lammps.dat ] 
then
	echo "Can't find input file. Are you sure you are in the right folder?"
	exit 1
fi

python ../../src/lammps-to-oxDNA.py init_lammps.dat
diff_lines=$(cat <(diff $CORRECT_OUTPUT $OUTPUT) <(diff $CORRECT_TOP $OUTPUT_TOP))

if [ $? -ne 0 ]
then
	echo "TEST FAILED";
	exit 1
else
	echo "TEST PASSED";
	#rm $OUTPUT
	exit 0
fi 
