#!/bin/bash

CORRECT_OUTPUT_CONF="correct_open_output.oxdna"
CORRECT_OUTPUT_TOP="correct_open_output.top"
OUTPUT_CONF="centerline_open.dat.oxdna"
OUTPUT_TOP="centerline_open.dat.top"

if [ ! -s centerline_open.dat ]
then
	echo "Can't find input files. Are you sure you are in the right folder?"
	exit 1
fi

rm $OUTPUT_CONF $OUTPUT_TOP 2> /dev/null
python ../../src/XYZ_oxDNA.py centerline_open.dat -p 0.1 --dsDNA --open --seed=1000
(diff $CORRECT_OUTPUT_CONF $OUTPUT_CONF > /dev/null) && (diff $CORRECT_OUTPUT_TOP $OUTPUT_TOP > /dev/null)

if [ $? -ne 0 ]
then
	echo "TEST FAILED";
	exit 1
else
	echo "TEST PASSED";
	rm $OUTPUT_CONF $OUTPUT_TOP
	exit 0
fi 
