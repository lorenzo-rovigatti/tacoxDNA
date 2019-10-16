#!/bin/bash

CORRECT_OUTPUT="correct_output.dat"
CORRECT_TOP="correct_output.top"
OUTPUT_CONF="input.ma.oxdna"
OUTPUT_TOP="input.ma.top"

if [ ! -s input.ma ] 
then
	echo "Can't find the input file. Are you sure you are in the right folder?"
	exit 1
fi

rm $OUTPUT_CONF $OUTPUT_TOP 2> /dev/null
python ../../src/vHelix_oxDNA.py input.ma
(diff $CORRECT_OUTPUT $OUTPUT_CONF > /dev/null) && (diff $CORRECT_TOP $OUTPUT_TOP > /dev/null)

if [ $? -ne 0 ]
then
	echo "TEST FAILED";
	exit 1
else
	echo "TEST PASSED";
	rm $OUTPUT_CONF $OUTPUT_TOP
	exit 0
fi 
