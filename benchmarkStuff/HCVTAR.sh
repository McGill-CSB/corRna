#!/bin/bash


echo "Running benchmark script for HCV sequence"
./serverbench.sh AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCCGACCCCGC 
cp  b-results.txt HCV-results.txt

echo "Running benchmark script for TAR sequence"
./serverbench.sh GGUCUCUCUGGUUAGACCAGAUCUGAGCCUGGGAGCUCUCUGGCUAACUAGAGAACC 
cp  b-results.txt TAR-results.txt



#echo "Running benchmark script for HCV sequence"
#./benchmark.sh AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCCGACCCCGC 
#cp  b-results.txt HCV-results.txt

#echo "Running benchmark script for TAR sequence"
#./benchmark.sh GGUCUCUCUGGUUAGACCAGAUCUGAGCCUGGGAGCUCUCUGGCUAACUAGAGAACC 
#cp  b-results.txt TAR-results.txt
