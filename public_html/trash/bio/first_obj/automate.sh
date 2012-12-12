#!/bin/bash
#for (( i=0; i <= 1000; i++))
#do
#generate the RNA...
	rm bpp-results.txt
	cat hepac.dat @- | RNAfold -p -d3
#done
#parses all the files
	./ps-adapter.pl
##Prossess the files...

echo "Correlating against wildtype_dp.ps.in";

for f in *in
do 
	if [ "wildtype_dp.ps.in" != $f ];
	then
		python correlate.py -x wildtype_dp.ps.in -y $f
	fi
done

echo "Check bpp-results.txt for comparing against all base bair probability";
