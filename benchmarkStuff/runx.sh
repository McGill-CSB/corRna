#!/bin/bash

#first agrument is the sequence
#secon argument is the number of mutations

function runx {
	touch run.txt
	echo "Running $1 with $2 mutation(s)"
	echo "m=$2" > run.txt
	for i in `seq 1 $3`
	do
		echo "Running trial: $i"
		./script.py -b -r -n 10 -s $1 -m $2
	done
}

