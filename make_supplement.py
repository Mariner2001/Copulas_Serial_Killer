import csv, re, math
def pairs(res):   # T_i, D_i as printed by copula_check.py
    out = []
    for l in open(res, encoding='utf-8'):
        m = re.match(r'\s*(\d+)\s+(.+?)\s+->\s+(.+?)\s+(\d+)\s+([0-9.]+)\s*$', l)
        if m: out.append((int(m.group(4)), float(m.group(5))))
    return out
def fmt_date(d):
    y, m, dd = d.split('-'); return f"{dd}/{m}/{y}"
def esc(x): return x.replace('&', r'\&').replace('_', r'\_')
# ---------------------------------------------------------------- S1 Sutcliffe
ma = "A A M M A M M M A M A M M M M M M A A M".split()
names = ["Anna Rogulskyj","Olive Smelt","Wilma McCann","Emily Jackson","Marcella Claxton","Irene Richardson","Patricia Atkinson","Jayne MacDonald",
         "Maureen Long","Jean Jordan","Marilyn Moore","Yvonne Pearson","Helen Rytka","Vera Millward","Josephine Whitaker","Barbara Leach",
         "Marguerite Walls","Upadhya Bandara","Theresa Sykes","Jacqueline Hill"]
S = list(csv.DictReader(open('sutcliffe.csv'))); P = pairs('res_sutcliffe.txt')
s1 = []
for i, r in enumerate(S):
    T, D = (P[i] if i < len(P) else ('---', None))
    s1.append(f"{i+1} & {names[i]} & {fmt_date(r['date'])} & {ma[i]} & {float(r['lat']):.3f}, ${float(r['lon']):.3f}$ & {T} & {'---' if D is None else f'{D:.1f}'}\\\\")
# ---------------------------------------------------------------- S2 Chikatilo
imp = {"Dunenkova": "July 1983", "Kutsyuba": "July 1983", "Unknown woman 1983": "July--Sept.\\ 1983", "Chuchulina": "after 19/09/1983",
       "Tsana": "c.\\ 28/07/1984", "Unknown woman Tashkent": "8--11/08/1984", "Unknown woman 1988": "1--4/04/1988"}
C = list(csv.DictReader(open('chikatilo.csv'))); CS = list(csv.DictReader(open('sites_chikatilo.csv'))); P = pairs('res_chikatilo.txt')
s2 = []
for i, (r, c) in enumerate(zip(C, CS)):
    T, D = (P[i] if i < len(P) else ('---', None))
    mark = r"$^\dagger$" if r['label'] in imp else ""
    s2.append(f"{i+1} & {esc(r['label'])} & {fmt_date(r['date'])}{mark} & {esc(c['locality'])} & {float(r['lat']):.3f}, {float(r['lon']):.3f} & {T} & {'---' if D is None else f'{D:.1f}'}\\\\")
# ---------------------------------------------------------------- S3 Ridgway
appr = {"Sherrill": "20/10--07/11/1982", "Brockman": "about 24/12/1982", "Rois": "31/05--13/06/1983", "Barczak": "October 1986",
        "Hayes": "1987", "Reeves": "1990", "Yellow Robe": "1998"}
RR = list(csv.DictReader(open('ridgway_disposal.csv'))); RS = list(csv.DictReader(open('sites_ridgway.csv'))); P = pairs('res_ridgway.txt')
s3 = []
for i, (r, c) in enumerate(zip(RR, RS)):
    T, D = (P[i] if i < len(P) else ('---', None))
    mark = r"$^\dagger$" if r['label'] in appr else ""
    s3.append(f"{i+1} & {esc(r['label'])} & {fmt_date(r['date'])}{mark} & {esc(c['encounter_area'])} & {esc(c['disposal_site'])} & {float(r['lat']):.3f}, ${float(r['lon']):.3f}$ & {T} & {'---' if D is None else f'{D:.1f}'}\\\\")
# ---------------------------------------------------------------- S4 sensitivity
def tau_line(label, res):
    txt = open(res, encoding='utf-8').read()
    n = int(re.search(r'pairs n = (\d+)', txt).group(1)); tau = float(re.search(r'Kendall tau = (-?[0-9.]+)', txt).group(1))
    z = tau * math.sqrt(9*n*(n-1)/(2*(2*n+5))); p = math.erfc(abs(z)/math.sqrt(2))
    return f"{label} & {n} & ${tau:.3f}$ & ${p:.3f}$\\\\"
s4 = [tau_line("Sutcliffe, main analysis", "res_sutcliffe.txt"), tau_line("Chikatilo, main analysis", "res_chikatilo.txt"),
      tau_line("\\quad custody (13/09--12/12/1984) deducted", "sens_chikatilo_custody.txt"),
      tau_line("\\quad 43 convictions upheld in 1993", "sens_chikatilo_43upheld.txt"),
      tau_line("Ridgway, main analysis (body-recovery sites)", "res_ridgway.txt"),
      tau_line("\\quad encounter (last-sighting) areas", "sens_ridgway_encounter.txt"),
      tau_line("\\quad dates of other public lists", "sens_ridgway_disposal_wikidates.txt"),
      tau_line("\\quad including Rebecca Marrero (49th victim)", "sens_ridgway_disposal_marrero.txt")]
doc = r"""\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[english]{babel}
\usepackage{amsmath,amssymb}
\usepackage{booktabs,longtable}
\usepackage{xcolor}
\usepackage[hidelinks]{hyperref}
\usepackage{geometry}
\geometry{a4paper,top=2.2cm,bottom=2.2cm,left=2cm,right=2cm}
\renewcommand{\thetable}{S\arabic{table}}
\setlength{\parindent}{0pt}\setlength{\parskip}{0.5em}
\title{\textbf{Supplementary material}\\[0.3em]\large Copula-based analysis of space-time dependence in serial homicide: evidence from three high-volume, long-duration series}
\date{}
\begin{document}
\maketitle

\section*{S1. Data, sources and conventions}
Each series lists the offences in chronological order; $T_i$ (days) and $D_i$ (km) are the interval and the haversine distance ($R=6{,}371$~km) from offence $i$ to offence $i+1$. The files \texttt{sutcliffe.csv}, \texttt{chikatilo.csv} and \texttt{ridgway\_disposal.csv} (columns \texttt{label, date, lat, lon}) and the script \texttt{copula\_check.py} reproduce every table and figure of the paper; variants used in the sensitivity analyses (Table~S4) are also provided. The script \texttt{revision\_checks.py} reproduces Tables~S5 and~S6 and the bootstrap intervals of Table~7 of the paper.

\textbf{Sutcliffe (Table S1).} Attack dates and street-level coordinates of the 20 attacks, as in the original chronology \citep[from][]{cross,bilton}. Displacements are computed from the coordinates; the displacement column of the first version of the paper did not match them and has been corrected.

\textbf{Chikatilo (Table S2).} The 52 murders of which he was convicted (Laura Sarkisyan, whose murder he was cleared of at trial, is excluded); the double murder of Tatyana and Svetlana Petrosyan counts as one offence, giving 51 offences and 50 pairs. Dates and localities were compiled from the victim table of the English-language Wikipedia article on Chikatilo, which draws on Cullen (1993), Krivich and Ol'gin (1993) and Conradi (1992). Coordinates are those of the locality (town, settlement, park or railway station) and are approximate (errors of up to about 10~km). Dates marked $\dagger$ are imprecise in the sources and were fixed as follows, preserving the order of the sources: July 1983 (two victims) $\to$ 10 and 20 July; ``July--September 1983'' (unknown woman, between Gudkov and Chuchulina) $\to$ 1 September; ``after 19 September 1983'' $\to$ 19 September; ``c.\ 28 July 1984'' $\to$ 28 July; ``8--11 August 1984'' $\to$ 9 August; ``1--4 April 1988'' $\to$ 2 April. Other lists give 31 July 1985 for Pokhlistova and 8 March 1989 for Ryzhova. Five pairs of consecutive offences share a locality ($D_i=0$). The rank-based analyses use these zeros as tied smallest values; for the parametric fits they are set to 1~km, below the smallest positive displacement (5.1~km), and Table~S5 shows the effect of other values.

\textbf{Ridgway (Table S3).} The 44 victims with known dates among the 48 of the 2003 plea (the four unidentified victims are excluded). Dates of last sighting, last-sighting areas and body-recovery sites follow the list released at the plea (Associated Press, 5 November 2003), based on the prosecutor's summary (Maleng, 2003). Offences are located at the body-recovery site (coordinates of the recovery cluster, approximate); for Denise Bush, whose remains were found partly south of Seattle and partly in Tigard, Oregon, the first site is used. Last-sighting areas, known only at corridor level, are used in a sensitivity analysis. Dates marked $\dagger$ are approximate in the plea list and were fixed as follows: Sherrill (20 October--7 November 1982) $\to$ 20 October; Brockman (about 24 December 1982) $\to$ 24 December; Rois (31 May--13 June 1983) $\to$ 2 June; Barczak (October 1986) $\to$ 17 October; Hayes (1987, last seen leaving a Portland jail) $\to$ 7 February 1987; Reeves (1990) $\to$ 5 March 1990; Yellow Robe (1998) $\to$ 15 January 1998, the last four from other public lists. Five pairs of consecutive offences share a recovery cluster ($D_i=0$; four victims at the Green River and three at Star Lake). They are treated in the same way (smallest positive displacement 2.0~km; Table~S5).

\begin{table}[htbp]
\centering\small
\caption{Peter Sutcliffe: chronology of the 20 attacks (M\,=\,murder; A\,=\,attempted murder). Coordinates in decimal degrees; $T_i$ and $D_i$ lead from each attack to the next.}
\begin{tabular}{clcccrr}
\toprule
\# & Victim & Date & M/A & Latitude, Longitude & $T_i$ (d) & $D_i$ (km)\\
\midrule
""" + "\n".join(s1) + r"""
\bottomrule
\end{tabular}
\end{table}

{\small
\begin{longtable}{rllp{3.6cm}crr}
\caption{Andrei Chikatilo: chronology of the 51 offences (52 convictions). $\dagger$: imprecise date (see S1).}\\
\toprule
\# & Victim & Date & Locality & Latitude, Longitude & $T_i$ (d) & $D_i$ (km)\\
\midrule
\endfirsthead
\toprule
\# & Victim & Date & Locality & Latitude, Longitude & $T_i$ (d) & $D_i$ (km)\\
\midrule
\endhead
""" + "\n".join(s2) + r"""
\bottomrule
\end{longtable}}

{\footnotesize
\begin{longtable}{rllp{3.1cm}p{3.1cm}crr}
\caption{Gary Ridgway: chronology of the 44 dated victims of the 2003 plea. Coordinates are those of the body-recovery site. $\dagger$: approximate date (see S1).}\\
\toprule
\# & Victim & Last seen & Last-sighting area & Body-recovery site & Latitude, Longitude & $T_i$ (d) & $D_i$ (km)\\
\midrule
\endfirsthead
\toprule
\# & Victim & Last seen & Last-sighting area & Body-recovery site & Latitude, Longitude & $T_i$ (d) & $D_i$ (km)\\
\midrule
\endhead
""" + "\n".join(s3) + r"""
\bottomrule
\end{longtable}}

\begin{table}[htbp]
\centering\small
\caption{Kendall's $\hat\tau$ between intervals and displacements in the main and sensitivity analyses (asymptotic two-sided $p$-values).}
\begin{tabular}{lrrr}
\toprule
Series and variant & $n$ & $\hat\tau$ & $p$\\
\midrule
""" + "\n".join(s4) + r"""
\bottomrule
\end{tabular}
\end{table}

\section*{S2. Additional analyses (second revision)}
\textbf{Zero displacements (Table S5).} A recorded displacement of zero means that two consecutive offences were assigned the same locality or recovery cluster; the true distance is small but unknown. Rank-based analyses treat these zeros as tied smallest values and do not depend on any value assigned to them. The log-normal distribution gives no probability to zero, so for the parametric fits (Table~4 of the paper) and the IFM estimates that depend on them (Table~6 and Figure~1) the zeros were set to 1~km. Table~S5 repeats these fits with values from 0.1 to 5~km and with the zeros treated as left-censored. With every replacement value the log-normal remains preferred to the exponential and the IFM estimates change little, and the censored fits give similar log-normal parameters; no conclusion of the paper depends on this choice.

\textbf{Stability of the trigger (Table S6).} Bootstrap intervals for the operational parameters are given in Table~7 of the paper (2,000 resamples of the pairs $(T_i,D_i)$; seed 2026). Table~S6 shows how the empirical 85th percentile of the intervals, and three times the median interval, would have evolved as offences were linked. With 14 intervals (15 offences), the probability that all of them lie below the true 85th percentile is $0.85^{14}\approx0.10$; in the three series the estimates available after 15 offences differed from those of the complete series by factors of up to 2.6, whereas three times the median varied less.

""" + open('tables_S5_S6.tex', encoding='utf-8').read() + r"""

\clearpage
\section*{S3. References}
Bilton, M. (2003). \textit{Wicked Beyond Belief}. HarperCollins, London.\\
Conradi, P. (1992). \textit{The Red Ripper}. True Crime, London.\\
Cross, R. (1981). \textit{The Yorkshire Ripper}. Granada, London.\\
Cullen, R. (1993). \textit{The Killer Department}. Pantheon, New York.\\
Krivich, M. and Ol'gin, O. (1993). \textit{Comrade Chikatilo}. Barricade Books, Fort Lee.\\
Maleng, N. (2003). \textit{Prosecutor's Summary of the Evidence: State of Washington v.\ Gary Leon Ridgway}. King County Superior Court.\\
Associated Press (2003). Victims of Green River Killer, 5 November 2003.\\
Wikipedia (2026). Andrei Chikatilo --- Victims. \url{https://en.wikipedia.org/wiki/Andrei_Chikatilo} (accessed 24 September 2026).
\end{document}
"""
doc = doc.replace(r"\citep[from][]{cross,bilton}", "(from Cross, 1981, and Bilton, 2003)")
open('supplementary_material.tex', 'w', encoding='utf-8').write(doc)
print("rows:", len(s1), len(s2), len(s3))
