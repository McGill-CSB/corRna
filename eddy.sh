#!/bin/bash

touch b-results.txt
echo "Running $1"
echo "time, correlation, #sequences"
echo "m=1" > b-results.txt

echo "RNAMutants and FixedCGSampling" >> b-results.txt
echo "RNAMutants and FixedCGSampling" 
./script.py -r -g -n 20 -o -b -m 5 -s $1 > trash.txt
cat results.txt >> b-results.txt
echo "###########################" >> b-results.txt
echo "Structural Heuristic" >> b-results.txt
echo "Structural Heuristic" 
./backupsubopt.py -k -m 5 -n 20 -s $1 >> b-results.txt
echo "###########################" >> b-results.txt
echo "Mutation Heuristic" >> b-results.txt
echo "Mutation Heuristic" 
./wrapper.py -m 5 -n 20 -l 10 -t 1 -d -c -s $1 >> b-results.txt
