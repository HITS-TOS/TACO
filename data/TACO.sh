#!/bin/bash

echo "--------------------------"
echo "---- run TACO ----"
echo "--------------------------"

#MAXLWD=0.

for dir in ./0*.dir
do
    cd $dir
    #if [ -f peaksMLE.csv ]
    #then
	#echo "peaksMLE.csv already present. Skipping $dir"
    #else
        echo "$dir"
        echo "filter lightcurve"
        time Rscript ../../src/filter_lightcurve.R
        echo "compute pds"
        time python3 ../../src/pds.py
        echo "numax estimate"
        time Rscript ../../src/numax_estimate.R
        echo "background fit"
        time python3 ../../src/background_fit.py --bins 300
#       echo "background summary"
 #       python3 ../../../../lib/background_summary.py
        echo "find peaks"
        time Rscript ../../src/peakFind.R --snr 1.1 --prob 0.0001 --minAIC 2
        echo "MLE optimisation for peaks"
     	time Rscript ../../src/peaksMLE.R  --minAIC 2
        echo "mode ID 02"
        time Rscript ../../src/peakBagModeId02.R
        echo "find peaks"
        time Rscript ../../src/peakFind.R --snr 1.1 --prob 0.0001 --minAIC 2 --removel02 TRUE
        echo "MLE optimisation for peaks"
        time Rscript ../../src/peaksMLE.R  --minAIC 2 --removel02 TRUE --init peaksMLE.csv --mixedpeaks mixedpeaks.csv
        echo "MLE optimisation for peaks"
        time Rscript ../../src/peaksMLE.R  --minAIC 2 --finalfit TRUE --init peaksMLE.csv --mixedpeaks mixedpeaksMLE.csv

    #fi
    cd - > /dev/null
    #break
done
