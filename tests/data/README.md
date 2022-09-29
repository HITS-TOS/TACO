# Single execution of TACO modules

Rscript src/numax_estimate.R -p tests/data/test_numax_<n>_pds.csv -s tests/data/test_numax_<n>_summary.csv

python3 src/background_fit.py \
   --pds tests/data/test_background_fit/pds.csv \
   --summary tests/data/test_background_fit/summary.csv \
   --ofac_pds tests/data/test_background_fit/ofac_pds.csv \
   --bins 300 --seed 42

Rscript src/peakFind.R --snr 1.1 --prob 0.0001 --minAIC 2 --pds tests/data/test_peak_find_pds.csv --ofac_pds tests/data/test_peak_find_ofac_pds.csv --summary tests/data/test_peak_find_summary.csv

Rscript src/peaksMLE.R --minAIC 2 --pds tests/data/test_peaks_mle_pds_bgr.csv --init tests/data/test_peaks_mle_peaks.csv --summary tests/data/test_peaks_mle_summary.csv

cd tests/data/test_peak_bag_mode_id02
Rscript ../../../src/peakBagModeId02.R

cd tests/data/test_peak_bag_period_spacing
python3 ../../../src/peakBagPeriodSpacing.py
