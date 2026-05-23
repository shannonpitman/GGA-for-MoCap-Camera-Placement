# Cross-check: Conference paper vs OCP chapter

Compared `IEEE-conference-template-062824.tex` against `OCP.tex` and the
underlying simulation logs (`Results/7Cams/*.txt`).

Two classes of issue are reported separately: (1) **inconsistencies**
between the two documents (or between either document and the run
logs); (2) **content present in the chapter that should arguably be in
the paper**.

---

## 1. Inconsistencies

### 1.1 Mutation rate `μ` (high priority)

| Location | Value |
|---|---|
| Conference paper, Sec. V.E ("Selection, Mutation...") | `μ = 0.2` |
| Conference paper, Sec. VI (Experimental Setup) | `μ = 0.2` |
| Chapter, Table 6.4 (`tab:ga-params`) | `μ = 0.5` |
| Chapter prose ("$0.5 \times 42 \approx 21$ genes per offspring") | `μ = 0.5` |
| `Results/7Cams/*.txt` headers | `Mutation Rate: 0.50` |

The simulation logs are the ground truth. The compressed results table
I just wrote assumes `μ = 0.5`. **Fix:** change both occurrences in the
conference paper to `μ = 0.5`.

### 1.2 Number of independent seeds per modality (high priority)

| Location | Value |
|---|---|
| Conference paper, Sec. VI | "$n = 30$ independent random seeds" per modality |
| Chapter and underlying logs (7Cams) | 40 runs per (TT,GM,start) cell, i.e. 80 runs per (TT,GM) cell at $k=7$; 240 7-cam runs in total |

The conference paper undercounts the campaign. **Fix:** change to "$n=40$
seeds per start-strategy per cell" (or "$80$ runs per cell pooled across
cold- and warm-start") to match the data. If the paper does not want to
introduce warm-start in the methodology, change Sec. VI to "$n=40$ random
seeds per cell, cold-start only" and use only cold-start subset numbers
in the table.

### 1.3 Mutation step `σ` interpretation (medium priority)

The conference paper, Sec. V.E states `σ = 0.1` "($0.1\,\text{m}$ for
positions, $5.7^\circ$ for orientations)". The chapter, Sec.
`sec:param-selection` states the same `σ = 0.1` but claims this
perturbs "by $\pm 0.2$~m and $\pm 0.2$~rad or $\pm 11.5^\circ$ at one
standard deviation".

Both can't be right: with `σ = 0.1`,
$1\sigma = 0.1\,\text{rad} \approx 5.73^\circ$, not $11.5^\circ$.
The chapter's number ($\pm 11.5^\circ$) corresponds to $2\sigma$. **Fix:**
correct the chapter's prose to "$\sigma = 0.1$ corresponds to
$\pm 5.7^\circ$ at one standard deviation; over the $\pm 2\sigma$
operating range this is $\pm 11.5^\circ$". The conference paper's
phrasing is correct.

### 1.4 Cost-evaluation count (low priority)

Conference paper, Sec. III ("Resolution Uncertainty"): "Embedded in a
GA that performs $> 22{,}000$ cost evaluations per run".

Chapter, Sec. `sec:termination`: "$N_\text{pop}+G_\text{max}\!\cdot\!N_C
= 100 + 100\!\times\!150 = 15{,}100$ evaluations".

Compressed conference results section also uses $15{,}100$.

The $22{,}000$ figure does not match $p_C = 1.0$ and the stated
parameters. **Fix:** change "$> 22{,}000$" to "$> 15{,}000$" (or
$15{,}100$ exactly) in the conference paper. If the paper is using a
different upper-bound counting convention (e.g. counting all candidate
chromosomes including those replaced by elitism), say so explicitly,
otherwise the chapter's $15{,}100$ is correct.

### 1.5 Citation key for Chen \& Davis 2008 (low priority)

Conference paper, Sec. II.A: "refined with extensive empirical
validation by the same authors in 2008 `\cite{chen_camera_2008}`".

Chapter uses `\textcite{chen_occlusion_2008}` throughout.

**Fix:** make the key consistent across both documents; the chapter's
key (`chen_occlusion_2008`) appears to be the one that resolves
correctly in the bib.

### 1.6 Citation style: Corke RVC toolbox

Conference paper, Sec. VI cites `\cite{corkeRoboticsVisionControl2023}`.

Chapter, Sec. `sec:sim-environment` cites
`\parencite{corke_petercorkervc3-python_2026}` (year 2026 and a
distinct key).

The two cite different editions of the same toolbox. **Fix:** pick the
edition that was actually used in the simulation code and use one
consistent key across both documents.

### 1.7 Elitism delay $it_e$ (medium priority)

Conference paper, Sec. V.E describes elitism as merge-sort-truncate
with no mention of the $50$-generation delay.

Chapter, Sec. `sec:elitism` explicitly: "The mechanism activates only
after the elitism delay $it_e = 50$ generations: during the
exploratory phase $g < it_e$ the offspring population fully replaces
the parent population."

This is an important methodological detail (it is what produces the
characteristic flattening of the convergence median at $g\approx 50$
that my compressed results section refers to). **Fix:** add one
sentence to Sec. V.E of the paper, e.g. "Elitism activates only after
a delay of $it_e=50$ generations; before that, offspring fully replace
parents, allowing the search to spread before exploitation is
enforced."

### 1.8 Warm-starting (high priority — methodology gap)

Conference paper does not mention warm starting anywhere.

Chapter, Sec. `sec:initial-pop` and `sec:sim-convergence`: warm
starting (seeding the initial population with elite chromosomes from
prior cold-start runs) is one of the two main empirical claims of the
campaign, reducing median cost by 19--30% and tightening run-to-run
variability by 44--72%.

My compressed conference results section mentions warm starting
briefly. **Fix one of two ways:**

- (preferred) Add a short paragraph to Sec. V or Sec. VI of the paper
  describing the warm-start mechanism. The simplest framing is: "An
  optional warm-start mode replaces a small fraction of the
  perimeter-initialised population with elite chromosomes from prior
  cold-start runs." Then keep the warm-start sentence in the results.
- (alternative) Restrict the conference paper's results to cold-start
  runs only and remove the warm-start sentence from the compressed
  results.

### 1.9 Lens types in headline run (low priority — terminology)

Conference paper, Sec. VI: "four narrow-angle cameras ... and three
wide-angle cameras".

Chapter, Sec. `sec:chromosome-rep`: "three wide-angle cameras and four
narrow-angle cameras", but then "cameras at even-indexed slots are
assigned wide-angle lenses, yielding three wide-angle and four
narrow-angle cameras".

These agree numerically. Both consistent with the OptiTrack rig in
Duncan McMillan. No change required.

---

## 2. Pertinent content present in the chapter but missing from the paper

These are not strictly inconsistencies but they affect whether the
paper's claims are defensible without the chapter alongside.

### 2.1 OptiTrack hardware specification (high priority)

Chapter, Sec. `sec:hardware-specs` includes Table~\ref{tab:optitrack-specs}
with sensor resolution ($1280\times1024$), pixel pitch
($\rho = 4.57\,\mu\text{m}$, back-solved from datasheet), and FOV
($56^\circ$ narrow / $\sim 82^\circ$ wide). The pixel-pitch back-solve
is consequential: the chapter notes that passing the default
$\rho = 10\,\mu\text{m}$ from Corke's `CentralCamera` class would
simulate a $\sim 99^\circ$ horizontal FOV and silently inflate every
camera's coverage.

**Recommendation:** add a short subsection or 2-line note to Sec. VI
of the paper stating the back-solved pixel pitch and FOV values used.
Without this, a reviewer cannot reproduce the cost numbers.

### 2.2 Baseline-angle distribution figure (medium priority)

Chapter, Sec. `sec:sim-costeval`, Fig.~\ref{fig:baseline-angles}
(`BaselineAngles_GAvsOpti_UAV.pdf`,
`BaselineAngles_GAvsOpti_UGV.pdf`). The "in-band fraction" headline
numbers ($73.2\%$ vs $64.8\%$ etc.) are quoted in the compressed
conference results but not visually supported. If column space allows,
including this figure (it is the cleanest single-figure summary of the
geometric mechanism) would meaningfully strengthen Sec. III of the
paper's results.

### 2.3 Cost-field decomposition figure (medium priority)

Chapter, Sec. `sec:sim-costeval`, Fig.~\ref{fig:costfield}
(`CostField_GAvsOpti_*.pdf`). Shows per-point uncertainty and
occlusion across the volume. Important because it demonstrates that
the GA's cost reduction is spatial-uniform, not just an averaging
effect. If column space allows, include in compressed form.

### 2.4 Grid-spacing sensitivity figure (medium priority)

Chapter, Sec. `sec:sim-sensitivity`, Fig.~\ref{fig:spacing-sensitivity-uav}.
The conference results section asserts "$\Delta=1\,\text{m}$ deviates
by $\leq 5\%$ from the $\Delta=0.1\,\text{m}$ reference whilst
delivering a $\sim 250\times$ speedup". This is a methodological
defence that a reviewer is likely to query. Including either the
figure or a small table would head off the question.

### 2.5 Cross-modality evaluation (low priority)

Chapter, Sec. `sec:sim-optconfig` reports that the UAV-best chromosome
evaluated on the UGV grid yields $\approx 0.034$ and vice-versa
$\approx 0.072$, both well below the ad-hoc baseline but $\sim 2\times$
above the modality-matched best.

The conference paper currently claims "a single algorithm produces a
UAV-tuned ... and a UGV-tuned ... configuration by altering only the
target-point grid" as a contribution (Sec. I, item 4) but never shows
that the two are meaningfully different. One sentence in the results
section would close this loop.

### 2.6 Population diversity figure (low priority)

Chapter, Sec. `sec:sim-convergence` uses
Fig.~\ref{fig:diversity-7c} as the diagnostic that distinguishes
healthy elitist convergence from premature convergence. The
compressed conference results refers to this argument in prose
("diversity collapses by generation $\sim 50$ whilst the best cost
continues to drop") but cannot show it. Optional in the paper; not
load-bearing if column space is tight.

### 2.7 Per-target coverage table (low priority)

Chapter, Sec. `sec:sim-costeval`, Table~\ref{tab:coverage-stats}.
The compressed conference results quotes the headline numbers in
prose; if a second table fits in Sec. VII.D of the paper, including
it would let the reader verify the $0\%$-unobserved claim without
needing the chapter.

---

## 3. Summary

| Severity | Issue | Action |
|---|---|---|
| High | `μ` mismatch (paper: 0.2, chapter+logs: 0.5) | Update paper to 0.5 |
| High | Seed count mismatch ($n=30$ in paper vs $40$/$80$ in chapter+logs) | Update paper |
| High | Warm starting not mentioned in paper methodology | Add 1--2 sentences to Sec. V/VI |
| High | Pixel pitch back-solve / FOV not reported in paper | Add to Sec. VI |
| Medium | `σ`-to-degree conversion typo in chapter ($11.5^\circ$ → $5.7^\circ$) | Update chapter |
| Medium | Cost-evaluation count ($22{,}000$ vs $15{,}100$) | Update paper |
| Medium | Elitism delay $it_e$ omitted from paper | Add 1 sentence |
| Medium | Baseline-angle / cost-field figures missing from paper | Optional include |
| Medium | Grid-spacing sensitivity justification missing from paper | Optional include |
| Low | Chen \& Davis 2008 citation key inconsistent | Pick one key |
| Low | Corke RVC toolbox key inconsistent | Pick one key |
| Low | Cross-modality claim not demonstrated in paper | Add 1 sentence to results |
