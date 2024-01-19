#!/bin/bash

CORRECT_OUTPUT="correct_output.dat"
OUTPUT="ds.dat.lammps"

if [ ! -s ds.dat ]  || [ ! -s ds.top ]
then
	echo "Can't find input files. Are you sure you are in the right folder?"
	exit 1
fi

rm $OUTPUT 2> /dev/null
python3 ../../src/tacoxDNA/oxDNA_LAMMPS.py ds.top ds.dat
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
