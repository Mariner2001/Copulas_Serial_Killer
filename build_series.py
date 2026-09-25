"""Builds the chronologies of the three series (one CSV per series: label,date,lat,lon).
Coordinates are approximate (locality/cluster level) and must be checked against the primary sources."""
import csv
# ---------------------------------------------------------------- Chikatilo (52 convictions; Sarkisyan, acquitted, excluded)
C = {  # locality-level coordinates (approx.)
 "Shakhty": (47.7085, 40.2160), "Kirpichnaya st.": (47.735, 40.160), "Donleskhoz st.": (47.840, 40.100),
 "Novoshakhtinsk": (47.7575, 39.9364), "Novocherkassk": (47.411, 40.104), "Kazachi Lagerya": (47.500, 40.050),
 "Donskoy (Bagayevsky)": (47.390, 40.280), "Bagayevskaya": (47.320, 40.390), "Rostov (centre)": (47.222, 39.719),
 "Rostov (Aviators' Park)": (47.258, 39.790), "Rostov (Botanical Garden)": (47.234, 39.656), "Rostov (Don bank)": (47.210, 39.735),
 "near Rostov (Aksay)": (47.270, 39.870), "Krasny Sulin": (47.883, 40.067), "Krasnodar airport": (45.034, 39.170),
 "Enem (Adygea)": (44.925, 38.905), "Tashkent": (41.300, 69.240), "Domodedovo (Vostryakovo)": (55.410, 37.760),
 "Revda": (56.798, 59.908), "Zaporizhzhia": (47.838, 35.140), "Leningrad": (59.939, 30.316), "Ilovaisk": (47.925, 38.201),
 "Kolchugino": (56.300, 39.383)}
chik = [  # (label, date, locality)  dates: Wikipedia (EN) victim table; imprecise dates -> documented convention
 ("Zakotnova", "1978-12-22", "Shakhty"), ("Tkachenko", "1981-09-03", "Rostov (Don bank)"), ("Biryuk", "1982-06-12", "Donskoy (Bagayevsky)"),
 ("Volobuyeva", "1982-07-25", "Krasnodar airport"), ("Pozhidayev", "1982-08-13", "Enem (Adygea)"), ("Kuprina", "1982-08-16", "Kazachi Lagerya"),
 ("Karabelnikova", "1982-09-08", "Shakhty"), ("Kuzmin", "1982-09-15", "Shakhty"), ("Stalmachenok", "1982-12-11", "Novoshakhtinsk"),
 ("Dunenkova", "1983-07-10", "Rostov (Aviators' Park)"), ("Kutsyuba", "1983-07-20", "Shakhty"), ("Gudkov", "1983-08-09", "Rostov (Aviators' Park)"),
 ("Unknown woman 1983", "1983-09-01", "Novoshakhtinsk"), ("Chuchulina", "1983-09-19", "Kirpichnaya st."), ("Shevkun", "1983-10-27", "Shakhty"),
 ("Markov", "1983-12-27", "Novocherkassk"), ("Shalapinina", "1984-01-09", "Rostov (Aviators' Park)"), ("Ryabenko", "1984-02-21", "Rostov (Aviators' Park)"),
 ("Ptashnikov", "1984-03-24", "Novoshakhtinsk"), ("T. and S. Petrosyan", "1984-05-25", "Kirpichnaya st."), ("Bakulina", "1984-06-22", "Bagayevskaya"),
 ("Illarionov", "1984-07-10", "Rostov (centre)"), ("Lemesheva", "1984-07-19", "Kirpichnaya st."), ("Tsana", "1984-07-28", "Rostov (Aviators' Park)"),
 ("Golosovskaya", "1984-08-02", "Rostov (Aviators' Park)"), ("Alekseyeva", "1984-08-07", "Rostov (Don bank)"), ("Unknown woman Tashkent", "1984-08-09", "Tashkent"),
 ("Seydaliyeva", "1984-08-13", "Tashkent"), ("Chepel", "1984-08-28", "Rostov (Don bank)"), ("Luchinskaya", "1984-09-06", "Rostov (Aviators' Park)"),
 ("Pokhlistova", "1985-08-01", "Domodedovo (Vostryakovo)"), ("Gulyayeva", "1985-08-27", "Shakhty"), ("Makarenkov", "1987-05-16", "Revda"),
 ("Bilovetsky", "1987-07-29", "Zaporizhzhia"), ("Tereshonok", "1987-09-15", "Leningrad"), ("Unknown woman 1988", "1988-04-02", "Krasny Sulin"),
 ("Voronko", "1988-05-15", "Ilovaisk"), ("Muratov", "1988-07-14", "Donleskhoz st."), ("Ryzhova", "1989-02-28", "Shakhty"),
 ("Dyakonov", "1989-05-11", "Rostov (centre)"), ("Moiseyev", "1989-06-20", "Kolchugino"), ("Varga", "1989-08-19", "near Rostov (Aksay)"),
 ("Khobotov", "1989-08-28", "Shakhty"), ("Kravchenko", "1990-01-14", "Shakhty"), ("Makarov", "1990-03-07", "Rostov (Botanical Garden)"),
 ("Zuyeva", "1990-04-04", "Donleskhoz st."), ("Petrov", "1990-07-28", "Rostov (Botanical Garden)"), ("Fomin", "1990-08-14", "Novocherkassk"),
 ("Gromov", "1990-10-17", "Donleskhoz st."), ("Tishchenko", "1990-10-30", "Kirpichnaya st."), ("Korostik", "1990-11-06", "Donleskhoz st.")]
# ---------------------------------------------------------------- Ridgway (44 dated victims of the 2003 plea; 4 undated Jane Does excluded)
R = {  # body-recovery sites (cluster level, approx.)
 "Green River, Kent": (47.390, -122.250), "S of Sea-Tac airport": (47.422, -122.300), "N of Sea-Tac airport": (47.470, -122.310),
 "S Seattle suburbs": (47.470, -122.340), "Federal Way": (47.322, -122.312), "N Seattle": (47.720, -122.335), "Tukwila": (47.475, -122.260),
 "Tigard, OR": (45.431, -122.772), "Sumner": (47.203, -122.240), "Star Lake": (47.357, -122.289), "Auburn cemetery": (47.293, -122.195),
 "Auburn": (47.305, -122.215), "Maple Valley": (47.393, -122.047), "Enumclaw": (47.204, -121.992), "Lake Sawyer": (47.335, -122.036),
 "SR18 / I-90": (47.500, -121.820), "Snoqualmie": (47.520, -121.800), "SR410, 12 mi E of Enumclaw": (47.170, -121.750),
 "SR410 E of Enumclaw": (47.180, -121.850), "SR410 near Enumclaw": (47.195, -121.920), "Kent (Cottonwood Park)": (47.405, -122.260),
 "I-90 near North Bend": (47.470, -121.720), "I-90 E of North Bend": (47.430, -121.640), "Seward Park, Seattle": (47.550, -122.253),
 "SR18 near GRCC, Auburn": (47.312, -122.180), "Des Moines Way S": (47.450, -122.305)}
E = {  # last-sighting (encounter) areas (approx.)
 "Pacific Hwy S (Sea-Tac strip)": (47.443, -122.296), "Tacoma": (47.253, -122.444), "Aurora Ave N": (47.720, -122.345),
 "downtown Seattle": (47.607, -122.335), "International District": (47.598, -122.324), "south Seattle": (47.530, -122.290),
 "south Seattle suburbs": (47.480, -122.310), "Northgate": (47.706, -122.325), "Portland, OR (jail)": (45.515, -122.677)}
P, T_, AV, DS, SS, SE, NG, PO, ID, TA = ("Pacific Hwy S (Sea-Tac strip)", "Tacoma", "Aurora Ave N", "downtown Seattle", "south Seattle suburbs",
                                          "south Seattle", "Northgate", "Portland, OR (jail)", "International District", "Tacoma")
ridg = [  # (label, date [AP 2003 plea list; Wikipedia where AP gives only month/year], disposal site, encounter area)
 ("Coffield", "1982-07-08", "Green River, Kent", T_), ("Lovvorn", "1982-07-17", "S of Sea-Tac airport", P), ("Bonner", "1982-07-25", "Green River, Kent", P),
 ("Chapman", "1982-08-01", "Green River, Kent", P), ("Hinds", "1982-08-11", "Green River, Kent", P), ("Mills", "1982-08-12", "Green River, Kent", P),
 ("Milligan", "1982-08-29", "Star Lake", P), ("Meehan", "1982-09-15", "S Seattle suburbs", P), ("Estes", "1982-09-20", "Federal Way", P),
 ("Rule", "1982-09-26", "N Seattle", AV), ("Bush", "1982-10-08", "Tukwila", P), ("Summers", "1982-10-09", "N of Sea-Tac airport", DS),
 ("Sherrill", "1982-10-20", "Tigard, OR", ID), ("Brockman", "1982-12-24", "Sumner", P), ("A. Smith", "1983-03-03", "Star Lake", P),
 ("Williams", "1983-03-08", "Star Lake", P), ("Mathews", "1983-04-10", "Star Lake", P), ("Childers", "1983-04-16", "S of Sea-Tac airport", P),
 ("Gabbert", "1983-04-17", "Star Lake", P), ("Pitsor", "1983-04-17", "Auburn cemetery", DS), ("Malvar", "1983-04-30", "Auburn", P),
 ("Christensen", "1983-05-03", "Maple Valley", P), ("Authorlee", "1983-05-22", "Enumclaw", P), ("Wims", "1983-05-23", "N of Sea-Tac airport", DS),
 ("Antosh", "1983-05-31", "Lake Sawyer", P), ("Rois", "1983-06-02", "Star Lake", SS), ("Naon", "1983-06-08", "S of Sea-Tac airport", P),
 ("Ware", "1983-07-19", "S Seattle suburbs", DS), ("Thompson", "1983-07-25", "SR18 / I-90", P), ("Buttram", "1983-08-18", "Snoqualmie", SS),
 ("Abernathy", "1983-09-05", "SR410, 12 mi E of Enumclaw", DS), ("Winston", "1983-09-12", "Kent (Cottonwood Park)", NG),
 ("Feeney", "1983-09-28", "I-90 near North Bend", DS), ("Bello", "1983-10-11", "SR410 E of Enumclaw", P), ("Avent", "1983-10-26", "SR410 near Enumclaw", DS),
 ("Plager", "1983-10-30", "I-90 E of North Bend", SS), ("Nelson", "1983-11-01", "I-90 near North Bend", P), ("Yates", "1983-12-23", "I-90 E of North Bend", SE),
 ("West", "1984-02-06", "Seward Park, Seattle", SE), ("C. Smith", "1984-03-21", "SR18 near GRCC, Auburn", P), ("Barczak", "1986-10-17", "Auburn", P),
 ("Hayes", "1987-02-07", "SR410 E of Enumclaw", PO), ("Reeves", "1990-03-05", "SR410 near Enumclaw", P), ("Yellow Robe", "1998-01-15", "Des Moines Way S", P)]
def write(fn, rows):
    with open(fn, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["label", "date", "lat", "lon"]); w.writerows(rows)
write("chikatilo.csv", [(l, d, *C[s]) for l, d, s in chik])
write("chikatilo_custody.csv", [(l, ("1985-05-03" if l == "Pokhlistova" else d), *C[s]) for l, d, s in chik])  # placeholder, see below
up = {"Zakotnova","Tkachenko","Pozhidayev","Stalmachenok","Shalapinina","Tsana","Bilovetsky","Voronko","Kravchenko"}
write("chikatilo_43upheld.csv", [(l, d, *C[s]) for l, d, s in chik if l not in up])
write("ridgway_disposal.csv", [(l, d, *R[s]) for l, d, s, e in ridg])
write("ridgway_encounter.csv", [(l, d, *E[e]) for l, d, s, e in ridg])
wiki = {"Williams": "1983-03-14", "Mathews": "1983-04-22", "Childers": "1983-04-14", "Pitsor": "1983-04-16", "Christensen": "1983-05-04",
        "Ware": "1983-07-18", "Buttram": "1983-08-23", "Sherrill": "1982-10-22"}
rw = sorted([(l, wiki.get(l, d), *R[s]) for l, d, s, e in ridg], key=lambda r: r[1])
write("ridgway_disposal_wikidates.csv", rw)
m = [(l, d, *R[s]) for l, d, s, e in ridg]; m.insert(13, ("Marrero", "1982-12-03", *R["Auburn"]))
write("ridgway_disposal_marrero.csv", m)
# supplementary listing (site names)
with open("sites_chikatilo.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["label", "date", "locality"]); w.writerows(chik)
with open("sites_ridgway.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["label", "date", "disposal_site", "encounter_area"]); w.writerows(ridg)
print(len(chik), "Chikatilo offences;", len(ridg), "Ridgway offences")
