#!/bin/bash

CORRECT_OUTPUT_CONF="correct_simple_output.oxdna"
CORRECT_OUTPUT_TOP="correct_simple_output.top"
OUTPUT_CONF="centerline.dat.oxdna"
OUTPUT_TOP="centerline.dat.top"
CONF_DIFF_BIN="python3 ../conf_diff.py"

if [ ! -s centerline.dat ] || [ ! -s sequence.dat ]
then
	echo "Can't find input files. Are you sure you are in the right folder?"
	exit 1
fi

rm $OUTPUT_CONF $OUTPUT_TOP 2> /dev/null
python3 ../../src/tacoxDNA/XYZ_oxDNA.py centerline.dat -p 0.05 -q sequence.dat --dsDNA --closed
($CONF_DIFF_BIN $CORRECT_OUTPUT_CONF $OUTPUT_CONF > /dev/null) && (diff $CORRECT_OUTPUT_TOP $OUTPUT_TOP > /dev/null)

if [ $? -ne 0 ]
then
	echo "TEST FAILED";
	exit 1
else
	echo "TEST PASSED";
	rm $OUTPUT_CONF $OUTPUT_TOP
	exit 0
fi
