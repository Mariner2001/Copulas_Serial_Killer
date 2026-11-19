DATA AND CODE (Supplementary Material)
"Copula-based analysis of space-time dependence in serial homicide: evidence from three
high-volume, long-duration series"

Main series (label,date,lat,lon):   sutcliffe.csv, chikatilo.csv, ridgway_disposal.csv
Location / site names:               sites_chikatilo.csv, sites_ridgway.csv
Sensitivity variants:                chikatilo_custody.csv (excluding the 90 days in custody in 1984)
                                     chikatilo_43upheld.csv (43 convictions upheld in 1993)
                                     ridgway_encounter.csv (last-sighting areas)
                                     ridgway_disposal_wikidates.csv (dates from other public lists)
                                     ridgway_disposal_marrero.csv (including the 49th victim)
Results (B = 1000):                  res_sutcliffe.txt, res_chikatilo.txt, res_ridgway.txt
Sensitivity results (B = 100):       sens_*.txt
Second-revision analyses:            revision_checks.txt, tables_S5_S6.tex

To reproduce (Python 3 with numpy and scipy; run in this folder):
  python3 copula_check.py chikatilo.csv --B 1000 --seed 2026 --zero_km 1.0
                                  (Tables 2-6 and Figure 1 of the paper, one series at a time)
  python3 operational.py          (operational parameters, Table 7)
  python3 revision_checks.py      (bootstrap intervals of Table 7 and Tables S5-S6; about 3 minutes)
  python3 build_series.py         (rebuilds the CSV files from the lists)
  python3 make_supplement.py      (writes supplementary_material.tex with Tables S1-S6;
                                   run revision_checks.py first)

Coordinates are approximate: locality level for Chikatilo and recovery cluster for Ridgway
(see Section S1 of the Supplementary Material).
