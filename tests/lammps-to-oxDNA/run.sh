#!/bin/bash

CORRECT_OUTPUT="correct_output.dat"
OUTPUT="init_lammps.dat.out"

if [ ! -s init_lammps.dat ] 
then
	echo "Can't find input file. Are you sure you are in the right folder?"
	exit 1
fi

python ../../src/lammps-to-oxDNA.py init_lammps.dat
diff_lines=$(diff $CORRECT_OUTPUT $OUTPUT)

if [ $? -ne 0 ]
then
	echo "TEST FAILED";
	exit 1
else
	echo "TEST PASSED";
	rm $OUTPUT
	exit 0
fi 
