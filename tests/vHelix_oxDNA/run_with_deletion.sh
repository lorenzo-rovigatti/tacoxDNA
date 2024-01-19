#!/bin/bash

INPUT="input_with_deletion.ma"
CORRECT_OUTPUT="correct_output_with_deletion.dat"
CORRECT_TOP="correct_output_with_deletion.top"
OUTPUT_CONF="input_with_deletion.ma.oxdna"
OUTPUT_TOP="input_with_deletion.ma.top"
CONF_DIFF_BIN="python3 ../conf_diff.py"

if [ ! -s $INPUT ] 
then
	echo "Can't find the input file. Are you sure you are in the right folder?"
	exit 1
fi

rm $OUTPUT_CONF $OUTPUT_TOP 2> /dev/null
python3 ../../src/tacoxDNA/vHelix_oxDNA.py -e 12345 $INPUT
# we don't check the topology in this case
($CONF_DIFF_BIN $CORRECT_OUTPUT $OUTPUT_CONF > /dev/null)

if [ $? -ne 0 ]
then
	echo "TEST FAILED";
	exit 1
else
	echo "TEST PASSED";
	rm $OUTPUT_CONF $OUTPUT_TOP
	exit 0
fi