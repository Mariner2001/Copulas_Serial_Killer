DATA AND CODE (Supplementary Material)

Main series (label,date,lat,lon):            sutcliffe.csv, chikatilo.csv, ridgway_disposal.csv
Location / site names:                        sites_chikatilo.csv, sites_ridgway.csv
Sensitivity variants:                         chikatilo_custody.csv (excluding the 90 days in prison from 1984),
chikatilo_43upheld.csv (43 convictions upheld in 1993),
ridgway_encounter.csv (last sighting locations),
ridgway_disposal_wikidates.csv (dates from other public lists),
ridgway_disposal_marrero.csv (including victim 49)
Results (B = 1000):                           res_sutcliffe.txt, res_chikatilo.txt, res_ridgway.txt
Sensitivity results (B = 100):                sens_*.txt

To reproduce / Run:
python3 copula_check.py chikatilo.csv --B 1000 --seed 2026 --zero_km 1.0
python3 operational.py          (operational table, Table 7)
python3 build_series.py         (rebuilds the CSVs from the lists)
python3 make_supplement.py      (generates Tables S1-S4)

The coordinates for Chikatilo (location) and Ridgway (recovery location, by clusters) are approximate and should be cross-checked against primary sources prior to submission (see Section S1 of the Supplementary Material).
Requirements: Python 3, numpy, scipy.