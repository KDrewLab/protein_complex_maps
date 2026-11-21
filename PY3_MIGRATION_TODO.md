# Python 3 / Google Colab migration — TODO

Scope: everything under `protein_complex_maps/` **except `plant_map_website/`** (handled on another branch).
Baseline (Python 3.8): 30 of 77 in-scope `.py` files fail to compile; 54 in-scope files have 2to3 changes.

Test coverage tags on each item:
- **[TESTED]** — exercised by an existing unit test in `protein_complex_maps/test/`
- **[TESTED*]** — exercised only if an external dependency/tool/network is available (see Section 0)
- **[UNTESTED]** — no existing test reaches this code

Unit-test status. Python 2 baseline = conda env `py2` (2.7.18, pandas 0.24, scipy 1.2, sklearn 0.20, Biopython 1.76); Python 3 = anaconda base (3.8.5).
External tools: ClusterONE `/home/kdrew/programs/cluster_one-1.0.jar`, MCL `/home/kdrew/local/bin/mcl` (Java 11).
**In scope: ClusterONE and MCL only.** CFinder, CliXO and agglomod (networkx/agglomcluster) clustering paths are out of scope for now — not installed, not ported, not tested.

| Test module | Python 2 (baseline) | Python 3 (today) | Notes |
|---|---|---|---|
| `test_complex_comparison` | ✅ 18/18 | ✅ 18/18 | prints, `map()` and `xrange` fixed; stable across hash seeds 0-4 |
| `test_split_complexes` | ✅ 3/3 | ✅ 3/3 | prints converted; `len(...)/2` → `//` (stable across hash seeds 0-4) |
| `test_merge_complexes` | ✅ 1/1 | ✅ 1/1 | prints converted |
| `test_features` | ✅ 2/2 | ✅ 2/2 | — |
| `test_clustering_parameter_optimization` | ✅ `testClustering` (ClusterONE→MCL); `testClustering2` out of scope | ✅ `testClustering` | Text-mode temp files + `universal_newlines=True`; Py3 output identical to Py2 baseline `[['A','B','C'],['D','E','F']]`. `testClustering2` is CFinder→agglomod (out of scope) |
| `test_shared_bait_feature` | ✅ 9/9 | ✅ 9/9 | rpy2 removed (FDR now statsmodels); order-dependent pair lookups fixed (stable across hash seeds 0-9); `scipy.misc.comb` → `scipy.special.comb` |
| `test_protein_util` | ❌ import fails | ❌ import fails | **Code bug**: `protein_util.py` is half-converted (`import urllib.request` breaks Py2; `Bio.Alphabet` breaks Py3). No baseline exists |

Python 2 reference outputs to match after conversion:
- All asserted values in the 24 passing tests above.
- `testClustering` (ClusterONE → MCL) returns `[['A', 'B', 'C'], ['D', 'E', 'F']]` (the test has no assertion; see Section 6).

How to run the clustering test without networkx/agglomcluster installed: a wrapper that inserts empty stub modules into `sys.modules` before running unittest (only the `agglomod` path uses them). Per decision, networkx/agglomcluster are **not** being installed.

---

## 0. Get the existing test suite runnable (do first)

- [ ] Add a `requirements-test.txt` / test extra so tests can be run in Colab (`numpy pandas scipy scikit-learn biopython requests mpmath statsmodels`)
- [ ] `util/protein_util.py:12` — remove `import Bio.Alphabet as ba` and any `ba.` usages (Biopython ≥1.78 removed it) — **[TESTED]** (blocks `test_protein_util` import)
- [ ] `test/test_protein_util.py` — `get_from_uniprot` uses the legacy `www.uniprot.org/uniprot/?format=tab` API, retired in 2022. Port to `rest.uniprot.org` and/or mark the test as a network test (`@unittest.skipUnless(...)`) — **[TESTED*]**
- [ ] `util/protein_util.py` — already partially Py3 (`urllib.request`, `print_function`); finish converting by hand rather than re-running 2to3. Currently imports under neither Python
- [ ] (deferred — not installing networkx/agglomcluster) `clustering/clustering_parameter_optimization.py:21-22` — optionally move `import networkx` / `import agglomcluster.agglomod` inside the `agglomod` branch (~line 583) so the module imports without them
- [x] `test/test_clustering_parameter_optimization.py` — tool paths now read from env vars `CLUSTERONE_JAR`, `MCL_BIN`, `CFINDER_EXE`, `CFINDER_LICENSE`, `CLUSTER_TEMP_DIR`, defaulting to the local install paths
- [ ] (out of scope — CFinder) `test/test_clustering_parameter_optimization.py` — `testClustering2` could `skipUnless` CFinder exists so the suite is green; currently errors with `OSError: No such file or directory`
- [x] Replaced rpy2 BH FDR correction with `statsmodels.stats.multitest.multipletests(pvals, method='fdr_bh')[1]`; rpy2 imports removed. Verified identical to R `p.adjust(method='BH')` (max abs diff 1.1e-16 on 555 p-values). statsmodels is available in both `py2` (0.9) and base py3 (0.12)
  - [x] `features/shared_bait_feature.py`
  - [x] `features/shared_bait_feature_wsampling.py`
  - [ ] **Pre-existing bug**: `shared_bait_feature.py` computes `pval_corr` / `neg_ln_pval_corr` but the final column selection (`output_df[['gene_id1','gene_id2','neg_ln_pval','pair_count','pval']]`) drops them, so `bh_correct=True` has no effect on output. Decide whether to return them
  - [ ] **Pre-existing bug**: in `shared_bait_feature_wsampling.py` the `if bh_correct:` block is after `return output_dict` (unreachable)
  - [ ] `testSharedBaitFeature_bhcorrect` asserts only uncorrected values (use `get_pair()`) — add assertions on corrected values once the column is returned
- [ ] (deferred — not currently needed) `postprocessing_util/diffusion_clustering.py` — only remaining rpy2 user; embeds R code (igraph walktrap / hclust dendrogram, dendextend, dplyr). **Decide: keep R via rpy2 (Colab ships R + rpy2) or rewrite in Python** (python-igraph walktrap + scipy `cut_tree`; results won't be identical). If keeping R: — **[UNTESTED]**
  - [ ] `from itertools import izip` → `zip`
  - [ ] rpy2 3.x: `pandas2ri.activate()` deprecated, `pandas2ri.ri2py_dataframe` removed → use `with localconverter(default_converter + pandas2ri.converter): ...` and `rpy2py`
  - [ ] Move CRAN package install (runs at import time, lines 15-22) into `main()` or a setup step
- [ ] Run the suite with `python -m unittest discover protein_complex_maps/test` in a fresh Colab-like env and record baseline

## 1. Mechanical syntax fixes (2to3 / pyupgrade)

- [ ] Run `2to3 -w -n` (or `python-modernize`) on everything outside `plant_map_website/`; review the diff
- [ ] `print x` / `print >> f, x` → `print(...)` (~44 files)
  - [x] `evaluation/complex_comparison.py` (8), `preprocessing_util/complexes/complex_merge.py` (4), `preprocessing_util/complexes/split_complexes.py` (5) — all single-argument, so Py2-compatible without `__future__`; commented-out prints left as-is
  - **[TESTED]** in `complex_comparison.py`, `complex_merge.py`, `split_complexes.py`, `clustering_parameter_optimization.py`
  - **[UNTESTED]** everywhere else
- [ ] Implicit relative imports → absolute `protein_complex_maps.…` or `from . import …` (notably `complex_map_website/` → `complex_db`) — **[UNTESTED]**
- [ ] `xrange` → `range` (incl. `model_fitting/cross_validation/train_leaveout_divide.py`) — **[UNTESTED]**
- [ ] `except E, e:` → `except E as e:` (2 sites) — **[UNTESTED]**
- [ ] `.iteritems()` → `.items()` (`features/ExtractFeatures/canned_scripts/print_elution_info.py:11`) — **[UNTESTED]**
- [ ] `.next()` → `next(...)`, `itertools.izip` → `zip`, `cPickle` → `pickle` — **[UNTESTED]**
- [ ] Update shebangs `#!/usr/bin/python` / `#!/usr/bin/env python` → `python3` (incl. test files and `features/ExtractFeatures/*`)
- [ ] Remove now-redundant `from __future__ import print_function` (optional, cosmetic)

## 2. Semantic changes (review by hand)

### Dict views and lazy map/filter/zip
- [ ] Review every `.keys()/.values()/.items()` that is indexed, `len()`'d, iterated twice, or mutated during iteration (~51 lines / 18 files) — 2to3 wraps in `list()` conservatively
- [ ] Review every `map/filter/zip` result that is indexed or reused (~38 lines / 15 files)
  - **[TESTED]** for occurrences in `complex_comparison.py` (31 2to3 changes), `split_complexes.py`, `complex_merge.py`, `clustering_parameter_optimization.py`
  - **[UNTESTED]** elsewhere

- [x] **Silent bug** `evaluation/complex_comparison.py:377, 394, 469, 475, 517` — `np.max(map(...))` (always truthy on Py3) → `any(...)`. Same results; short-circuits (faster). Only difference: empty input returns False instead of raising — **[TESTED]**
- [x] `evaluation/complex_comparison.py:194, 202, 291, 292` — `map(len, ...)` → list comprehensions (line 202 would also have silently divided by 0 as the iterator was consumed twice) — **[TESTED]**
- [x] `evaluation/complex_comparison.py` — all `xrange` → `range` (446, 503 **[TESTED]**; 572, 579 in `chunks`/`rand_combinations` **[UNTESTED]**)

### Integer division (`/` now returns float)
- [x] `preprocessing_util/complexes/split_complexes.py:96` `split_id = len(complexes)/2` → `//` (slice index → TypeError) — **[TESTED]** (`testSplitComplexes*`)
- [x] `preprocessing_util/complexes/split_complexes.py:269` same, in `remove_intersection_random` — **[UNTESTED]**
- [ ] `util/protein_util.py:66` `range(len(ids)/ACC_QUERY_LENGTH + 1)` → `//` — **[TESTED*]** (network)
- [ ] `util/protein_util.py:110` same, in `get_from_uniprot_by_genename` — **[UNTESTED]**
- [ ] `evaluation/complex_comparison.py:578` `len(l)/n` in `rand_combinations` (function appears unused — consider deleting) — **[UNTESTED]**
- [ ] `features/ExtractFeatures/Features.py:212` `self.df + (1/len(self.df))` — Py2 added 0, Py3 adds a real pseudocount. **Decide intended behaviour — results will change** — **[UNTESTED]** (test only covers `sum_difference`, `pearsonR`)
- [ ] Grep remaining `/` between ints in non-`1.0*` expressions and confirm intent

### Bytes vs str
- [x] `clustering/clustering_parameter_optimization.py` `cluster_helper` ClusterONE + MCL path — **[TESTED]** by `testClustering`, Py3 output identical to Py2:
  - input temp files (lines 447, 528) `NamedTemporaryFile(mode='w')` (default binary mode rejects str on Py3)
  - ClusterONE `Popen(..., universal_newlines=True)` so `communicate()` returns str (`text=True` needs 3.7+; this also works on Py2)
  - MCL output file read with `"r"` instead of `"rb"` (MCL stdout isn't parsed, so its `Popen` is unchanged)
- [ ] `clustering/clustering_parameter_optimization.py:119` `main()` opens gold standard with `"rb"` → `"r"` — **[UNTESTED]**
- [ ] `clustering/clustering_parameter_optimization.py:448` (pre-existing leak) `dirTemp = tf.mkdtemp()` is created on every `cluster_helper` call but never removed; add `shutil.rmtree(dirTemp)` to the `finally` (respecting `nodelete`) — **[UNTESTED]**
- [ ] (out of scope — CFinder/CliXO) `clustering/clustering_parameter_optimization.py:479-507` CFinder `open(..., 'rb')` (487) and CliXO `communicate()` (498, 503) — same fix — **[UNTESTED]**
- [ ] `clustering/network_clustering.py:258` (ClusterONE) and MCL ~332 — same fix — **[UNTESTED]** (line 291 is CliXO — out of scope)
- [ ] `postprocessing_util/complex_enrichment.py:44,48` — CSV opened `'rb'` → `'r', newline=''` — **[UNTESTED]**
- [ ] `postprocessing_util/complex_enrichment.py:103` — drop `.encode('utf-8')` (yields `b'...'` in output) — **[UNTESTED]**
- [ ] Any Py2-produced `.pkl` inputs: load with `pickle.load(f, encoding='latin1')` (`complex_enrichment.py:54` checkpoint) — **[UNTESTED]**

## 3. Removed / changed library APIs

- [x] `scipy.misc.comb` → `scipy.special.comb` (`import scipy.misc as misc` → `import scipy.special as special`; 13 call sites). Verified: `test_shared_bait_feature` 9/9 on Py2 + Py3; `test_complex_comparison` 18/18 and `testClustering` pass on Py2 (Py3 still blocked by print statements)
  - [x] `evaluation/complex_comparison.py` (`clique_comparison`, `clique_comparison_exact`) — **[TESTED]**
  - [x] `features/shared_bait_feature.py` (`pval_old`) — **[TESTED]**
  - [x] `features/shared_bait_feature_wsampling.py` — **[UNTESTED]** (compiles on Py3)
- [ ] pandas `.ix` → `.iloc`
  - `evaluation/complex_comparison.py:165,177,193,211` (only on `topN != None` branch) — **[UNTESTED]** (no test passes `topN`)
  - `model_fitting/cross_validation/train_leaveout_divide.py:36-39` — **[UNTESTED]**
- [ ] `sklearn.metrics.ranking._binary_clf_curve` → `sklearn.metrics._ranking` or reimplement (`evaluation/plots/prcurve.py:17`) — **[UNTESTED]**
- [ ] `import matplotlib.mlab` — remove or replace usages in `evaluation/plots/prcurve.py`, `precision_score_curve.py`, `plot_chromatogram.py`, `plot_sparklines.py`, `model_fitting/cross_validation/evaluate_crossvalidation.py` — **[UNTESTED]**
- [ ] `pd.read_table` → `pd.read_csv(sep=...)` in `features/add_label.py`, `alphabetize_pairs.py`, `combine_feature_matrices.py` (deprecation, not breaking) — **[UNTESTED]**
- [ ] Sweep for other pandas ≥2 removals (`DataFrame.append`, `.iteritems`, `sort(columns=)`) after the run — **[UNTESTED]**

## 4. Packaging & Colab usability

- [ ] Add `pyproject.toml` (installable via `pip install git+https://github.com/...`) with core deps and optional extras: `r` (rpy2 — only `diffusion_clustering.py`, not currently needed), `agglomod` (networkx, agglomcluster — out of scope for now), `web` (flask, sqlalchemy, wtforms), `plots` (matplotlib, seaborn, bokeh)
- [ ] Add missing `__init__.py`: `postprocessing_util/`, `model_fitting/cross_validation/`, `model_fitting/svm_utils/`, `features/ExtractFeatures/canned_scripts/` (and `model_fitting/` itself)
- [ ] Replace hardcoded paths with CLI args / env vars:
  `clustering/network_clustering.py`, `clustering/clustering_parameter_optimization.py`, `evaluation/psweep_comparison.py`, `util/cluster2node_table_cdm_mod.py`, `util/label_training_interactions.py`, `util/label_interaction_scores.py`, `test/test_clustering_parameter_optimization.py`
- [ ] Document installing external tools in Colab (ClusterONE jar via Java, MCL via `apt-get install mcl` or build) — e.g. an example notebook cell
- [ ] Add an example Colab notebook covering the main pipeline
- [ ] Update `README.md` with Python 3 install/usage instructions

## 5. `complex_map_website/` (decide scope)

- [ ] Decide: port now, defer, or drop. It's a Flask app (not useful on Colab); 5 of the 30 non-compiling files live here
- [ ] If porting: print fixes, relative `complex_db` imports, check `requirements.txt` versions (Flask/SQLAlchemy APIs) — **[UNTESTED]**

## 6. Close the test gaps (recommended, small)

Only 7 of 77 in-scope modules have any test; none of the CLI scripts are tested. Cheapest high-value additions:
- [x] **Fixed flaky `test_shared_bait_feature` tests**: added `get_pair()` helper that matches a pair in either orientation and asserts exactly one row; all 25 active pair lookups use it. Previously Py2 always failed 3 tests and Py3 failed on some hash seeds. (Underlying output orientation is still hash-order dependent — optionally sort `gene_id1`/`gene_id2` in `shared_bait_feature` output for reproducible files)
- [ ] `testClustering` — add `assertEqual(sorted(map(sorted, predicted_clusters)), [['A','B','C'],['D','E','F']])` (Py2 baseline output); currently it only prints
- [ ] `complex_comparison` — add a test passing `topN=` (covers the `.ix` fixes)
- [ ] `Features.extract_features(..., resampling="poisson_noise")` — pin down the pseudocount behaviour decision
- [ ] `split_complexes.remove_intersection_random`
- [ ] `network_clustering` / `cluster_helper` with a stubbed subprocess (covers bytes/str without needing Java tools)
- [ ] `complex_enrichment` read/write round-trip on a tiny CSV
- [ ] `prcurve` smoke test on a tiny score file (covers `_binary_clf_curve`, `mlab`)
- [ ] Smoke test: `python -m py_compile` on every in-scope file (catches all Section 1 regressions cheaply)
- [ ] Optional: run a few key scripts end-to-end under Py2 and Py3 on the same small input and diff outputs
