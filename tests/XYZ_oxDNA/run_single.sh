#!/bin/bash

CORRECT_OUTPUT_CONF="correct_single_output.oxdna"
CORRECT_OUTPUT_TOP="correct_single_output.top"
OUTPUT_CONF="centerline.dat.oxdna"
OUTPUT_TOP="centerline.dat.top"

if [ ! -s centerline.dat ]
then
	echo "Can't find input files. Are you sure you are in the right folder?"
	exit 1
fi

rm $OUTPUT_CONF $OUTPUT_TOP 2> /dev/null
python ../../src/XYZ_oxDNA.py centerline.dat -p 0.1 -q sequence.dat --ssDNA --open
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
