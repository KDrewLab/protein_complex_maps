# Python 3 / Google Colab migration — TODO

Scope: everything under `protein_complex_maps/` **except `plant_map_website/`** (handled on another branch).
Baseline (Python 3.8): 30 of 77 in-scope `.py` files failed to compile at start; **now 0** (after 2to3 pass, Section 1).

Test coverage tags on each item:
- **[TESTED]** — exercised by an existing unit test in `protein_complex_maps/test/`
- **[TESTED*]** — exercised only if an external dependency/tool/network is available (see Section 0)
- **[UNTESTED]** — no existing test reaches this code

Unit-test status. Python 2 baseline = conda env `py2` (2.7.18, pandas 0.24, scipy 1.2, sklearn 0.20, Biopython 1.76); Python 3 = anaconda base (3.8.5).
External tools: ClusterONE `/home/kdrew/programs/cluster_one-1.0.jar`, MCL `/home/kdrew/local/bin/mcl` (Java 11).
**In scope: ClusterONE and MCL only.** CFinder, CliXO and agglomod (networkx/agglomcluster) clustering paths are out of scope for now — not installed, not ported, not tested.

| Test module | Python 2 (baseline) | Python 3 (today) | Notes |
|---|---|---|---|
| `test_complex_comparison` | ✅ 19/19 | ✅ 19/19 | prints, `map()` and `xrange` fixed; stable across hash seeds 0-4; `testComplexComparisonTopN` added |
| `test_split_complexes` | ✅ 3/3 | ✅ 3/3 | prints converted; `len(...)/2` → `//` (stable across hash seeds 0-4) |
| `test_merge_complexes` | ✅ 1/1 | ✅ 1/1 | prints converted |
| `test_features` | ✅ 2/2 | ✅ 2/2 | — |
| `test_clustering_parameter_optimization` | ✅ `testClustering` (ClusterONE→MCL); `testClustering2` out of scope | ✅ `testClustering` | Text-mode temp files + `universal_newlines=True`; Py3 output identical to Py2 baseline `[['A','B','C'],['D','E','F']]`. `testClustering2` is CFinder→agglomod (out of scope) |
| `test_shared_bait_feature` | ✅ 9/9 | ✅ 9/9 | rpy2 removed (FDR now statsmodels); order-dependent pair lookups fixed (stable across hash seeds 0-9); `scipy.misc.comb` → `scipy.special.comb` |
| `test_protein_util` | ❌ import fails (Py3-only module) | ✅ 1/1 (network) | Ported `get_from_uniprot` to UniProt REST API; skips if rest.uniprot.org unreachable. Module is Py3-only (no Py2 baseline ever existed) |

Python 2 reference outputs to match after conversion:
- All asserted values in the 24 passing tests above.
- `testClustering` (ClusterONE → MCL) returns `[['A', 'B', 'C'], ['D', 'E', 'F']]` (the test has no assertion; see Section 6).

How to run the clustering test without networkx/agglomcluster installed: a wrapper that inserts empty stub modules into `sys.modules` before running unittest (only the `agglomod` path uses them). Per decision, networkx/agglomcluster are **not** being installed.

---

## 0. Get the existing test suite runnable (do first)

- [x] ~~Add a `requirements-test.txt` / test extra~~ — not needed: tests only use `unittest` plus the core dependencies, which `pip install` brings in via `pyproject.toml`
- [x] `util/protein_util.py` — removed `Bio.Alphabet`; `get_sequences_uniprot(seqrecord=True)` now sets `annotations={'molecule_type': 'protein'}` on the `SeqRecord` — **[TESTED]** (import) / checked manually (SeqRecord path)
- [x] `util/protein_util.py` `get_from_uniprot` ported to `https://rest.uniprot.org/uniprotkb/stream` (legacy API retired 2022) — **[TESTED]**:
  - `requests` instead of `urllib` (str not bytes); query `accession:A OR accession:B`; `stream` endpoint so results aren't truncated at 25 (`search` paginates)
  - `UNIPROT_FIELD_MAP` translates legacy keywords callers use (`length`, `organism`, `genes`, `genes(PREFERRED)`, `protein+names`, `sequence`, `reviewed`, `id`) to REST field names — all checked manually against live API
  - only unique, well-formed accessions are queried (`UNIPROT_ACC_REGEX`); the REST API fails the **whole batch** with HTTP 400 if any accession is malformed, so malformed IDs are skipped with a warning. Isoforms (`P12345-2`) still return under the base accession, as before
  - checked 250 accessions across 3 batches: all returned
- [x] `test/test_protein_util.py` — `skipUnless` rest.uniprot.org is reachable (verified: skips cleanly with network blocked)
- [ ] `util/protein_util.py` `get_from_uniprot_by_genename` — still uses legacy URL; port to REST (`reviewed:yes` → `reviewed:true`, `organism:` → `organism_name:`/`organism_id:`) — **[UNTESTED]**, not used by in-scope modules
- [ ] `util/protein_util.py` `map_protein_ids` — legacy `/mapping/` endpoint is gone (HTTP 405). Port to async ID-mapping API (`POST /idmapping/run` → poll `/idmapping/status/{jobId}` → `/idmapping/uniprotkb/results/stream/{jobId}`); DB names changed (e.g. `ACC` → `UniProtKB_AC-ID`). Used by `convert_ids.py`, `postprocessing_util/cluster2node_table.py`, `util/cluster2node_table_cdm_mod.py` — **[UNTESTED]**; add a network test
- [x] `clustering/clustering_parameter_optimization.py` and `clustering/network_clustering.py`: `import networkx` / `import agglomcluster.agglomod` moved into the agglomod branch of `cluster_helper`, so both modules (and `clustering/trim_clusters2threshold.py`, which imports one of them) import without those packages (needed for Colab, which has networkx but not agglomcluster). Verified: full suite in py3colab (neither package installed) passes 36/37 with no stubs; `testClustering` passes on Py3.8 and Py2; agglomod path itself still untested (agglomcluster not installed anywhere) — **[UNTESTED]**
- [x] `test/test_clustering_parameter_optimization.py` — tool paths now read from env vars `CLUSTERONE_JAR`, `MCL_BIN`, `CFINDER_EXE`, `CFINDER_LICENSE`, `CLUSTER_TEMP_DIR`, defaulting to the local install paths
- [ ] (out of scope — CFinder) `test/test_clustering_parameter_optimization.py` — `testClustering2` could `skipUnless` CFinder exists so the suite is green; currently errors with `OSError: No such file or directory`
- [x] Replaced rpy2 BH FDR correction with `statsmodels.stats.multitest.multipletests(pvals, method='fdr_bh')[1]`; rpy2 imports removed. Verified identical to R `p.adjust(method='BH')` (max abs diff 1.1e-16 on 555 p-values). statsmodels is available in both `py2` (0.9) and base py3 (0.12)
  - [x] `features/shared_bait_feature.py`
  - [x] `features/shared_bait_feature_wsampling.py`
  - [x] **Pre-existing bug (fixed)**: `shared_bait_feature.py` computed `pval_corr` / `neg_ln_pval_corr` but the final column selection dropped them, so `--bh_correct` had no effect. Now BH is applied after the duplicate-pair removal (each pair is computed from both genes' sides, so 52 rows for 26 pairs on the test data) and the two columns are returned when `bh_correct=True`; output without `bh_correct` unchanged. (Exact duplication of every test turns out not to change BH values, but correcting the unique pairs is correct in general.) Verified: `pval_corr` matches R `p.adjust(method='BH')` exactly; uncorrected columns identical to the original — **[TESTED]**
  - [x] **Pre-existing bug (fixed)**: in `shared_bait_feature_wsampling.py` the `if bh_correct:` block was after `return output_dict` (unreachable); moved before it. Verified against R `p.adjust` exactly (each pair once there). `main()` calls it with the default `bh_correct=False`, so CLI output unchanged — **[UNTESTED]** by the suite
  - [x] `testSharedBaitFeature_bhcorrect` now asserts `pval_corr` / `neg_ln_pval_corr` for 4 pairs (expected values from R); `testSharedBaitFeature` asserts the corrected columns are absent without `bh_correct`
- [ ] (deferred — not currently needed) `postprocessing_util/diffusion_clustering.py` — only remaining rpy2 user; embeds R code (igraph walktrap / hclust dendrogram, dendextend, dplyr). **Decide: keep R via rpy2 (Colab ships R + rpy2) or rewrite in Python** (python-igraph walktrap + scipy `cut_tree`; results won't be identical). If keeping R: — **[UNTESTED]**
  - [ ] `from itertools import izip` → `zip`
  - [ ] rpy2 3.x: `pandas2ri.activate()` deprecated, `pandas2ri.ri2py_dataframe` removed → use `with localconverter(default_converter + pandas2ri.converter): ...` and `rpy2py`
  - [ ] Move CRAN package install (runs at import time, lines 15-22) into `main()` or a setup step
- [x] Run the suite in a fresh Colab-like env and record baseline — conda env `py3colab`: Python 3.12.11, numpy 2.5.5, pandas 3.1.1, scipy 1.18.1, scikit-learn 1.9.1, statsmodels 0.15.0, biopython 1.88, matplotlib 3.11.2, mpmath 1.3.0 (latest pip releases, 2026-11-21; may be newer than Colab's pins). Run with unittest discovery + networkx/agglomcluster stubs. **First run: 37 tests, 33 pass, 4 errors; after the fixes below: 36 pass, only `testClustering2` (CFinder, out of scope) errors:**
  - [x] `test_features` `testPearson`, `testSumDifference` — test `setUp` does `pd.DataFrame(d, dtype=float)` with a string `ids` column; old pandas silently skipped that column, new pandas raises `ValueError: could not convert string to float: 'a'`. Fixed in the test (build, `set_index('ids')`, then `astype(float)`); same DataFrame as before on pandas 1.1. Passes in py3colab, base Py3.8 and Py2
  - [x] `test_shared_bait_feature` `testSharedBaitFeatureBoundaryCaseOld` — `shared_bait_feature.py:53` `pval_old` uses `scipy.exp` (SciPy no longer re-exports numpy functions) → `np.exp` (fixed). Only such use in in-scope code. Passes in py3colab, base Py3.8 and Py2
  - `testClustering2` — CFinder, out of scope (same as before)
  - Static checks in `py3colab`: all 70 in-scope files compile with SyntaxWarnings as errors (Py3.12 invalid-escape warnings: none). `pd.read_table` still exists in pandas 3 (deprecation item below is cosmetic). `DataFrame.append` is gone: `features/shared_bait_feature_wsampling.py` `concat_df.append(tmp_df)` loop → one `pd.concat` (fixed). Verified by running the script end-to-end (numpy seeded) on a 60-row input: new code on Py3.8/pandas 1.1 gives byte-identical output and intermediate concat file to the original; new code on py3colab gives identical output; original fails on py3colab with `AttributeError` — **[UNTESTED]** by the suite (other `.append` hits are list appends). Chained-assignment grep (silently no-op under pandas 3 copy-on-write) found only dict-of-dict assignments

## 1. Mechanical syntax fixes (2to3 / pyupgrade)

- [x] Ran 2to3 on all 77 in-scope files (38 changed). **All 77 now compile on Python 3** (was 47/77). Method:
  - 50 files that already parsed as Py3 → `2to3 -w -n -p -x future` (`-p` so existing `print(a, b)` isn't turned into a tuple; `-x future` keeps `from __future__ import print_function`)
  - 27 Py2-only files → `2to3 -w -n`
  - reviewed dry-run diffs before writing
- [x] `print` statements → `print()` everywhere; 6 redundant `print((x))` from 2to3 tidied to `print(x)` (`evaluation/psweep_comparison.py`, `evaluation/plots/prcurve.py`)
- [x] `xrange` → `range`, `.iteritems()` → `.items()`, `except E, e` → `except E as e`, `cPickle` → `pickle`, `itertools.izip` → `zip`, `u''` prefixes dropped, `.next()` → `__next__`
- [x] **Real Py3 bugs fixed by 2to3's `list(...)` wrapping** — `df['col'] = map(...)` / `zip(...)` would store a single map object instead of a column: `features/combine_feature_matrices.py`, `features/shared_bait_feature.py`, `features/shared_bait_feature_wsampling.py`, `util/label_interaction_scores.py`, `postprocessing_util/pairwise2clusterid_multi.py`; `map(protid_set.add, ...)` used for side effects would silently do nothing: `postprocessing_util/cluster2node_table.py`, `util/cluster2node_table_cdm_mod.py`; `map(float, ...)` lists in `evaluation/psweep_comparison.py`
- [x] Tabs/spaces mix (TabError on Py3) in `clustering/network_clustering.py` and `preprocessing_util/ms1quant2elution.py` → leading tabs expanded at width 8 (Python 2's rule; verified identical Py2 AST before/after)
- [x] **Pre-existing bug** `model_fitting/svm_utils/svm_results2pairs.py:16` `except ValueError, TypeError:` caught only `ValueError` (binding it to the name `TypeError`) → `except (ValueError, TypeError):`
- [x] `evaluation/complex_comparison.py` `WeightedRandomGenerator`: 2to3 renamed `next` → `__next__`; added `next = __next__` alias so Py2 baseline tests still run — **[TESTED]**
- [ ] Noise from 2to3 (harmless, optional cleanup): `list(d.keys())` / `list(range(...))` wrapping in loops/comprehensions where a view/range is fine
- [ ] Implicit relative imports in `complex_map_website/` (`import complex_db`): 2to3 left them, works when run as a script from that directory; convert to absolute `protein_complex_maps.…` if imported as a package — **[UNTESTED]**
- [ ] Update shebangs `#!/usr/bin/python` / `#!/usr/bin/env python` → `python3` (incl. test files and `features/ExtractFeatures/*`)
- [ ] Remove now-redundant `from __future__ import print_function` (optional, cosmetic — kept for now so Py2 baseline comparisons still work)
- **Untested**: apart from the 7 tested modules, the 2to3 changes are only verified to compile, not to run

## 2. Semantic changes (review by hand)

### Dict views and lazy map/filter/zip
- [x] `.keys()/.values()/.items()` — 2to3 wrapped every non-loop use in `list()`
- [x] `map/filter/zip` — 2to3 wrapped every non-iteration use in `list()`; remaining unwrapped calls checked by hand: all are `multiprocessing` `Pool.map` (returns a list) or SQLAlchemy `.filter()` (unrelated)
  - **[TESTED]** in the 7 tested modules; **[UNTESTED]** (compile-only) elsewhere

- [x] **Silent bug** `evaluation/complex_comparison.py:377, 394, 469, 475, 517` — `np.max(map(...))` (always truthy on Py3) → `any(...)`. Same results; short-circuits (faster). Only difference: empty input returns False instead of raising — **[TESTED]**
- [x] `evaluation/complex_comparison.py:194, 202, 291, 292` — `map(len, ...)` → list comprehensions (line 202 would also have silently divided by 0 as the iterator was consumed twice) — **[TESTED]**
- [x] `evaluation/complex_comparison.py` — all `xrange` → `range` (446, 503 **[TESTED]**; 572, 579 in `chunks`/`rand_combinations` **[UNTESTED]**)

### Integer division (`/` now returns float)
- [x] `preprocessing_util/complexes/split_complexes.py:96` `split_id = len(complexes)/2` → `//` (slice index → TypeError) — **[TESTED]** (`testSplitComplexes*`)
- [x] `preprocessing_util/complexes/split_complexes.py:269` same, in `remove_intersection_random` — **[UNTESTED]**
- [x] `util/protein_util.py:66` `range(len(ids)/ACC_QUERY_LENGTH + 1)` → `//` — **[TESTED]**
- [x] `util/protein_util.py` same, in `get_from_uniprot_by_genename` — **[UNTESTED]**
- [x] `evaluation/complex_comparison.py:582` `len(l)/n` → `len(l)//n` in `rand_combinations` (unused in the repo; kept). Verified over 140 (len, n, m) cases with a fixed numpy seed: new code on Py3 yields exactly the same groups as the original on Py2 (original on Py3 differed in 37 cases, mostly fewer groups). Not covered by unit tests — **[UNTESTED]** by the suite
  - [ ] **Pre-existing quirk** (same on Py2/Py3): `chunks()` uses `i+n < len(l)`, so it also drops the last chunk when it is full (e.g. 4 items, n=2 → 1 chunk, not 2), and `rand_combinations` can yield fewer than `m` groups (it can also yield more — it never truncates to `m`). Should be `<=` if full chunks are intended. Warning comment added above `rand_combinations`
- [x] `features/ExtractFeatures/Features.py:214` `self.df + (1/len(self.df))` (poisson_noise pseudocount) — Py2 added 0 (integer division), Py3 would add a real 1/len(df). **Decision:** default reproduces Py2 (`1//len(df)`, i.e. no pseudocount) so results match previous runs; new `extract_features(..., real_pseudocount=True)` / `canned_scripts/extract_features.py --real_pseudocount` adds the real 1/len(df) (likely the original intent of the Traver-style code). Verified with seeded `--resampling poisson_noise` pearsonR runs on a 30x12 elution table: default on Py3 matches the original on Py2 (max diff 1.1e-16); `--real_pseudocount` on Py3 matches the original on Py3 exactly (the two settings differ by up to 0.35); default on Py2 unchanged. Unit tests don't cover resampling — **[UNTESTED]** by the suite
- [ ] Grep remaining `/` between ints in non-`1.0*` expressions and confirm intent

### Bytes vs str
- [x] `clustering/clustering_parameter_optimization.py` `cluster_helper` ClusterONE + MCL path — **[TESTED]** by `testClustering`, Py3 output identical to Py2:
  - input temp files (lines 447, 528) `NamedTemporaryFile(mode='w')` (default binary mode rejects str on Py3)
  - ClusterONE `Popen(..., universal_newlines=True)` so `communicate()` returns str (`text=True` needs 3.7+; this also works on Py2)
  - MCL output file read with `"r"` instead of `"rb"` (MCL stdout isn't parsed, so its `Popen` is unchanged)
- [x] `clustering/clustering_parameter_optimization.py:118` `main()` opens gold standard with `"r"` (done in the subprocess-bytes commit) — **[UNTESTED]**
- [x] (pre-existing leak) `cluster_helper` in `clustering/clustering_parameter_optimization.py` and `clustering/network_clustering.py`: `dirTemp = tf.mkdtemp()` was created on every call (only CFinder writes to it) but never removed. Added `shutil.rmtree(dirTemp, ignore_errors=True)` to the phase-1 `finally`, skipped with `nodelete`. Verified in an empty temp dir: `testClustering` (Py3 and Py2) and 3 `network_clustering` runs leave nothing (original left 1 dir per call); `nodelete=True` still keeps the files and dirs; clusters unchanged — **[TESTED]** by `testClustering` (runs the path; leftovers checked by hand)
- [ ] (out of scope — CFinder/CliXO) `clustering/clustering_parameter_optimization.py:479-507` CFinder `open(..., 'rb')` (487) and CliXO `communicate()` (498, 503) — same fix — **[UNTESTED]**
- [x] `clustering/network_clustering.py` `cluster_helper` ClusterONE + MCL path — same fix: input temp files (lines 232, 315) `NamedTemporaryFile(mode='w')`; ClusterONE `Popen(..., universal_newlines=True)`; MCL output was already read with `"r"`. Verified end-to-end (ClusterONE alone, ClusterONE+MCL, ClusterONE+MCL with `unweighted_one_step` + `trim2threshold`) on the unit-test network: Py3 output identical to original code on Py2; new code also unchanged on Py2; original code fails on Py3 with `TypeError: a bytes-like object is required`. No unit test — **[UNTESTED]** by the suite. Still open: CFinder `'rb'` (275) and CliXO `communicate()` (285) — out of scope; `dirTemp` leak fixed separately (see above)
- [x] `postprocessing_util/complex_enrichment.py:44,48` — opened `'r'` in the `"rb"` sweep; `newline=''` not needed (files are read with `readlines()`/`split()`, not the `csv` module)
- [x] `postprocessing_util/complex_enrichment.py:104` — `.encode('utf-8')` now only on Py2 (on Py3 it wrote `b'Caf\xc3\xa9...'` into the description column). Verified with a mocked g:Profiler response containing non-ASCII descriptions: Py3 output identical to original code on Py2 (after sorting columns — Py3 dicts keep insertion order, so the CSV column order differs, which is cosmetic); new code unchanged on Py2. No unit test — **[UNTESTED]** by the suite
- [x] **Text files opened with `"rb"`** → `"r"`: 45 reads in 23 files (only `"rb"` left is the pickle at `postprocessing_util/complex_enrichment.py:53`). `complex_map_website/load_db.py` (feeds `csv.reader`) uses `"r", newline=''`. Verified end-to-end: `prcurve.py` (pair IDs no longer print as `b'..'`; AUC/AP/precision-marker results identical to original Py2) and `complex_comparison.py` `main()` (deterministic metrics identical to original Py2; clique metrics are sampled and vary run-to-run in the original too). Unit tests unchanged — they don't reach these `main()` reads — mostly **[UNTESTED]**
- [x] **Text files opened for writing with `"wb"`** → `"w"`: 20 opens in 14 files (only `"wb"` left are the 2 pickles in `postprocessing_util/complex_enrichment.py`). All were plain `.write()` of str. Verified end-to-end vs original code on Py2 (outputs compared with lines/IDs order-normalised): `complex_merge.py`, `cluster2pairwise.py`, `cluster2cliquesize.py`, `trim_clusters2threshold.py` (run with networkx/agglomcluster stubbed) identical; `split_complexes.py` writes all 6 files correctly but see next item — **[UNTESTED]** by unit tests
- [x] **`split_complexes.py --random_seed` reproducible again, and identical on Py2 and Py3**: list order before each seeded `np.random.shuffle` came from `set` iteration, which Py3 randomises per process (numpy's seeded shuffle itself is identical across versions — checked numpy 1.16 vs 1.24). Set-derived lists now get a canonical order (`_canonical_order`) before shuffling: `merged_complexes` in `main()`, `large_pairs` via `sorted(set(clust))` (2 sites), and the intersection in `remove_intersection_random`. Caller-supplied lists and `complex_list` in `remove_overlapping_complexes` are already ordered and left alone; `remove_intersection` doesn't shuffle. Verified: identical splits on Py2 and Py3 (3 hash seeds) for default, `--size_threshold/--threshold_fraction` and `--merge_threshold 0.5`; `test_split_complexes`/`test_merge_complexes` pass on both under several hash seeds. **Note:** a given seed does not reproduce splits made with the old code (that order came from Py2's internal hash-table layout); keep old split files as the record — **[UNTESTED]** by unit tests (end-to-end checked)
- [x] Py2-produced `.pkl` inputs: `complex_enrichment.py` `--resume_from_checkpoint` now uses `load_checkpoint()` — tries a plain `pickle.load`, and on `UnicodeDecodeError` (Py2 pickle on Py3) reloads with `encoding='latin1'` and re-decodes the `description` column as UTF-8 (the Py2 script stored descriptions as UTF-8 bytes, which latin1 alone turns into mojibake `CafÃ©`). Verified with a checkpoint written by the original script on Py2 (mocked g:Profiler, non-ASCII descriptions): resume on Py3 gives output identical to resuming with the original code on Py2; Py3- and Py2-written checkpoints still resume correctly on their own version; original code on Py3 fails with `UnicodeDecodeError`. The only other pickle load is `features/ExtractFeatures/canned_scripts/average_features.py:27` `pd.read_pickle` — pandas already falls back to latin1 there (feature DataFrames, no pre-encoded strings), left as is; old-pandas pickles may still fail to unpickle on new pandas — **[UNTESTED]**. No unit test — **[UNTESTED]** by the suite

## 3. Removed / changed library APIs

- [x] `scipy.misc.comb` → `scipy.special.comb` (`import scipy.misc as misc` → `import scipy.special as special`; 13 call sites). Verified: `test_shared_bait_feature` 9/9 on Py2 + Py3; `test_complex_comparison` 18/18 and `testClustering` pass on Py2 (Py3 still blocked by print statements)
  - [x] `evaluation/complex_comparison.py` (`clique_comparison`, `clique_comparison_exact`) — **[TESTED]**
  - [x] `features/shared_bait_feature.py` (`pval_old`) — **[TESTED]**
  - [x] `features/shared_bait_feature_wsampling.py` — **[UNTESTED]** (compiles on Py3)
- [x] pandas `.ix` → `.loc` (label-based, matching what `.ix` did for these integer-labelled axes). Verified identical output: original `.ix` on Py2/pandas 0.24 vs `.loc` on Py3/pandas 1.1.3, for all four `topN` distributions (topN = None/1/2/4) and the CV train/leaveout splits
  - [x] `evaluation/complex_comparison.py:165,177,193,211` — **[UNTESTED]** by unit tests (no test or caller passes `topN`); verified by direct comparison
  - [x] `model_fitting/cross_validation/train_leaveout_divide.py:36-39` — **[UNTESTED]** by unit tests; verified by direct comparison
  - [x] **Pre-existing `topN` bugs (fixed)** `complex_comparison.py` — `topN` is never passed anywhere in the repo (checked all files, incl. positional calls), so no existing results change; `topN=None` paths unchanged:
    - off-by-one: label slice `.loc[:,:topN]` (as `.ix` did) kept topN+1 clusters → `.iloc[:,:topN]` (4 sites)
    - `mmr_pwmmr_hmean` accepted `topN` but ignored it → passes it through
    - `pwmmr` divided by all clusters (too small); `ppv` weighted by all clusters (NaN, so `acc` NaN); `precision_measure` used all cluster lengths (`ValueError`, so `precision_recall_product` too) → all now use the first `topN` clusters
    - **[TESTED]** by new `testComplexComparisonTopN`: for 9 metrics, `metric(topN=k)` equals the metric on a comparison built from `clusters[:k]` (2 datasets, every k incl. k > number of clusters), plus a hand-computed `mmr(topN=1)`; test fails on the previous code
- [x] `evaluation/plots/prcurve.py` `_binary_clf_curve`: `sklearn.metrics.ranking` removed in sklearn 0.24 → local copy of the function (adapted from sklearn, BSD-3) so it no longer depends on sklearn private internals. Verified identical `fps`/`tps`/`thresholds` vs sklearn 0.24 (`_ranking`) and 0.20 (`ranking`) on 300 random cases (ties, weights, single class, lists). End-to-end `prcurve.py` run on synthetic data (default and `--weight_positives`): AUC, average precision and every precision-marker recall/threshold identical to original code on Py2 — **[UNTESTED]** by unit tests
- [x] `import matplotlib.mlab` — was never actually broken (module still exists in matplotlib 3.x) and was unused in all 5 files; removed the dead imports from `evaluation/plots/prcurve.py`, `precision_score_curve.py`, `plot_chromatogram.py`, `plot_sparklines.py`, `model_fitting/cross_validation/evaluate_crossvalidation.py`
- [ ] `pd.read_table` → `pd.read_csv(sep=...)` in `features/add_label.py`, `alphabetize_pairs.py`, `combine_feature_matrices.py` (deprecation, not breaking) — **[UNTESTED]**
- [ ] Sweep for other pandas ≥2 removals (`DataFrame.append`, `.iteritems`, `sort(columns=)`) after the run — **[UNTESTED]**

## 4. Packaging & Colab usability

- [x] Add `pyproject.toml` (setuptools; name `protein_complex_maps`, version 0.1.0, MIT, `requires-python >=3.8`). Core deps: numpy, pandas, scipy, scikit-learn, statsmodels, biopython, mpmath, requests, matplotlib, seaborn (seaborn is core because `evaluation/plots/prcurve.py` etc. import it; bokeh is only used by the websites so not listed). Extras: `agglomod` (networkx; agglomcluster is not on PyPI), `r` (rpy2, python-igraph; `diffusion_clustering.py`). `plant_map_website` / `complex_map_website` excluded from the package. `.gitignore` now ignores `build/`, `dist/`, `*.egg-info/`. Verified: wheel contains all 16 in-scope packages (75 .py), no website files; `pip install <repo>` into a fresh venv (Python 3.12, deps from PyPI) then the suite run against the installed copy from outside the repo: 36/37 (only `testClustering2`, CFinder); `python -m protein_complex_maps.<module>` works for scripts; wheel also builds on Python 3.8
- [x] Add missing `__init__.py` (empty): `postprocessing_util/`, `model_fitting/`, `model_fitting/cross_validation/`, `model_fitting/svm_utils/`, `features/ExtractFeatures/canned_scripts/` — without them setuptools skipped those modules
- [x] Hardcoded tool paths: most hits in the listed files were commented-out lines. The real ones are the `clustering/network_clustering.py` and `clustering/clustering_parameter_optimization.py` CLI defaults. ClusterONE / MCL defaults now read `CLUSTERONE_JAR` / `MCL_BIN` (same names as the tests) and fall back to the previous values (ClusterONE fallback `/home/kdrew/programs/clusterone/cluster_one-1.0.jar` kept for now; MCL `/usr/bin/mcl` / `mcl`); explicit flags still win. CFinder defaults left (out of scope). Verified: resolved defaults with no env / env set / env + flags in both scripts; `network_clustering.py` run end-to-end with only the env vars gives the same 10 cluster files + params as the previous code given explicit flags
  - [x] **Pre-existing bug (fixed)**: both scripts did `os.path.exists(args.temp_dir)` before the `None` check, so the default `--temp_dir` (None) crashed with `TypeError` on Py2 and Py3; check order swapped
- [x] Document installing external tools in Colab — README section + notebook: `apt-get install mcl`; ClusterONE jar from https://paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar (checked: identical SHA-256 to the local jar; the notebook verifies the checksum)
- [x] Example notebook `notebooks/colab_quickstart.ipynb` (install, tool check, toy network → `network_clustering.py` ClusterONE→MCL → `ComplexComparison` metrics, optional test run). Executed locally end-to-end against the pip-installed package in a fresh Python 3.12 venv with the two `colab-setup` cells (pip/apt/wget, Colab paths) swapped for local paths; all cells run, tests 36/37 (only `testClustering2`). **Run on Google Colab 2026-11-23** (Python 3.12, OpenJDK 21, mcl 22-282): all cells run, same clusters and metrics as the local run (MMR 0.917, sensitivity 0.909, PPV 1.000), tests 36/37 (only `testClustering2`); harmless apt warning from Colab's own r2u (R packages) source. Covers clustering + evaluation only, not the earlier pipeline steps (elution profiles, features, model training)
- [x] `README.md`: "Python 3 / Google Colab" section (install, MCL/ClusterONE + `CLUSTERONE_JAR`/`MCL_BIN`, running scripts with `python -m`, tests, what is not ported). Test command checked from the repo root without `PYTHONPATH`: 36/37

## 5. `complex_map_website/` (decide scope)

- [ ] Decide: port now, defer, or drop. It's a Flask app (not useful on Colab); 5 of the 30 non-compiling files live here
- [ ] If porting: print fixes, relative `complex_db` imports, check `requirements.txt` versions (Flask/SQLAlchemy APIs) — **[UNTESTED]**

## 6. Close the test gaps (recommended, small)

Only 7 of 77 in-scope modules have any test; none of the CLI scripts are tested. Cheapest high-value additions:
- [x] **Fixed flaky `test_shared_bait_feature` tests**: added `get_pair()` helper that matches a pair in either orientation and asserts exactly one row; all 25 active pair lookups use it. Previously Py2 always failed 3 tests and Py3 failed on some hash seeds. (Underlying output orientation is still hash-order dependent — optionally sort `gene_id1`/`gene_id2` in `shared_bait_feature` output for reproducible files)
- [x] `testClustering` — now asserts the clusters (sorted) equal `[['A','B','C'],['D','E','F']]` (Py2 baseline output). Passes in py3colab, base Py3.8 and Py2
- [x] `complex_comparison` — `testComplexComparisonTopN` passes `topN=` (covers the `.ix` → `.iloc` fixes)
- [ ] `Features.extract_features(..., resampling="poisson_noise")` — test both `real_pseudocount` settings with a fixed seed
- [ ] `split_complexes.remove_intersection_random`
- [ ] `network_clustering` / `cluster_helper` with a stubbed subprocess (covers bytes/str without needing Java tools)
- [ ] `complex_enrichment` read/write round-trip on a tiny CSV
- [ ] `prcurve` smoke test on a tiny score file (covers `_binary_clf_curve`, `mlab`)
- [x] Smoke test: `test/test_compile.py` compiles every package `.py` file in memory (websites skipped; warnings as errors, so Py3.12 invalid escapes fail too). 75 files pass in py3colab, base Py3.8 and Py2; checked that it fails on a copy of the package with a planted `print "x"` file and a `"\d"` escape
- [ ] Optional: run a few key scripts end-to-end under Py2 and Py3 on the same small input and diff outputs
