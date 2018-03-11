#!/bin/bash

passed=0
tot=0
for f in $(ls -1 */*.sh)
do
	d=$(dirname $f)
	cd $d

	tot=$[tot + 1]
	
	bash $(basename $f) &> log.dat
	
	if [ $? -eq 0 ]
	then
		passed=$[passed + 1]
	else
		echo "$d: TEST FAILED"	
	fi
	
	cd ..
done

echo "$passed/$tot TESTS PASSED"
