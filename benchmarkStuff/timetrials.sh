#!/bin/bash

function runx {
        touch run.txt
        echo "Running $1 with $2 mutation(s)"
        echo "m=$2" > run.txt
        for i in `seq 1 $3`
        do
                echo "Running trial: $i"
                ./script.py -b -r -n 10 -s $1 -m $2 >> run.txt
        done
}

echo "\nn=60" >> results.txt
#Run n=60 for 1, 2, 3, 4, 5

runx UGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGG 4 4
cat run.txt >> results.txt

runx UGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGG 5 4
cat run.txt >> results.txt


#Run n=80 for 1, 2, 3, 4, 5

echo "\nn=80" >> results.txt
runx UGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGG 4 4
cat run.txt >> results.txt

runx UGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGG 5 4
cat run.txt >> results.txt


#Run n=100 for 1, 2, 3, 4, 5

echo "\nn=100" >> results.txt
runx UGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGG 4 4
cat run.txt >> results.txt

runx UGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGG 5 4
cat run.txt >> results.txt


#Run n=120 for 1, 2, 3, 4, 5

echo "\nn=120" >> results.txt
runx UGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGG 4 4
cat run.txt >> results.txt

runx UGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGGUGCCUGCCUCUUGGGAGGGG 5 4
cat run.txt >> results.txt

#echo "Running benchmark script for HCV sequence"
#./serverbench.sh AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCCGACCCCGC 
#cp  b-results.txt HCV-results.txt

#echo "Running benchmark script for TAR sequence"
#./serverbench.sh GGUCUCUCUGGUUAGACCAGAUCUGAGCCUGGGAGCUCUCUGGCUAACUAGAGAACC 
#cp  b-results.txt TAR-results.txt



#echo "Running benchmark script for HCV sequence"
#./benchmark.sh AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCCGACCCCGC 
#cp  b-results.txt HCV-results.txt

#echo "Running benchmark script for TAR sequence"
#./benchmark.sh GGUCUCUCUGGUUAGACCAGAUCUGAGCCUGGGAGCUCUCUGGCUAACUAGAGAACC 
#cp  b-results.txt TAR-results.txt
