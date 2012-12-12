#!/bin/bash

function runx {
        touch run.txt
        echo "Running $1 with $2 mutation(s)"
        echo "m=$2" > run.txt
        for i in `seq 1 $3`
        do
                echo "Running trial: $i"
                ./script.py -r -n 10 -s $1 -m $2 >> run.txt
        done
}

#Run n=20 for 1, 2, 3, 4, 5

runx UGCCUGCCUCUUGGGAGGGG 1 3
cat run.txt > results.txt


