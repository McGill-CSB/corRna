#!/bin/bash
touch b-results.txt

echo "Running $1 with 1 mutation"
echo "time, correlation, #sequences"
echo "m=1" > b-results.txt
for i in `seq 1 10`
do
	echo "run $i done"
	./script.py -b -r -n 10 -m 1 -s $1 >> b-results.txt
done


echo "Running $1 with 2 mutations"
echo "m=2" >> b-results.txt
for i in `seq 1 10`
do
	echo "run $i done"
	./script.py -b -r -n 10 -m 2 -s $1 >> b-results.txt
done


echo "Running $1 with 3 mutations"
echo "m=3" >> b-results.txt
for i in `seq 1 10`
do
	echo "run $i done"
	./script.py -b -r -n 10 -m 3 -s $1 >> b-results.txt
done


echo "Running $1 with 4 mutations"
echo "m=4" >> b-results.txt
for i in `seq 1 10`
do
	echo "run $i done"
	./script.py -b -r -n 10 -m 4 -s $1 >> b-results.txt
done


echo "Running $1 with 5 mutations"
echo "m=5" >> b-results.txt
for i in `seq 1 10`
do
	echo "run $i done"
	./script.py -b -r -n 10 -m 5 -s $1 >> b-results.txt
done
