protein_complex_maps
====================

Kdrew's scripts for handling protein complex map data

<!-- claude: section below written by claude for the py3-migration branch -->
Python 3 / Google Colab
-----------------------

The `py3-migration` branch runs on Python 3 (tested on 3.8 and 3.12 with current numpy/pandas/scipy). Migration status and remaining work: `PY3_MIGRATION_TODO.md`.

Example notebook: `notebooks/colab_quickstart.ipynb` [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/KDrewLab/protein_complex_maps/blob/py3-migration/notebooks/colab_quickstart.ipynb)

Install (Python >= 3.8):

    pip install git+https://github.com/KDrewLab/protein_complex_maps@py3-migration
    # or, from a clone of the repository:
    pip install .

External clustering tools (ClusterONE followed by MCL is the supported combination):

* MCL: `apt-get install mcl` (installs `/usr/bin/mcl`), or build from https://micans.org/mcl/
* ClusterONE: needs Java; download https://paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar
* Point the clustering scripts at them with environment variables, or the command-line options shown in brackets:

      export CLUSTERONE_JAR=/path/to/cluster_one-1.0.jar   # --clusterone_jar (network_clustering.py), --clusterone (clustering_parameter_optimization.py)
      export MCL_BIN=/usr/bin/mcl                          # --mcl_bin (network_clustering.py), --mcl (clustering_parameter_optimization.py)

Run scripts as modules, e.g.:

    python -m protein_complex_maps.clustering.network_clustering --input_network network.txt --output_file clusters.txt --twostep_combination clusterone mcl

Run the tests from the repository root (the clustering test uses `CLUSTERONE_JAR` / `MCL_BIN`; the UniProt test needs internet access; `testClustering2` needs CFinder and is expected to error):

    python -m unittest discover -s protein_complex_maps/test -t .

Not ported / not supported on this branch: CFinder, CliXO and agglomod clustering (agglomod needs `networkx` and `agglomcluster`), `postprocessing_util/diffusion_clustering.py` (R via rpy2), and the `complex_map_website` / `plant_map_website` Flask apps.

How to run:
1. Generate elution profiles  (preprocessing)
2. Calculate correlations (feature extraction)
3. Convert to pairwise files (feature extraction)
4. Build feature matrix (feature extraction)
4b. Ensure Common ID (feature extraction)
5. Split benchmark (Model Training/make_benchmark)
6. Train classifier (Model Training)
7. Predict interactions (Model Training)
8. Evaluate predicted interactions (Evaluation)
9. Cluster interactions (Clustering)
10. Evaluate predicted clusters (Evaluation)


Preprocessing  (src/preprocessing_util/)

Feature Extraction (src/features/)

Model Training (src/model_fitting/)
    SVM, LDA, tpot, other machine learning

Clustering (src/clustering/)

Evaluation (src/evalution/)


Filename conventions:

1. elution profiles .elut
2. tidy elution profiles .tidyelut
3. pairwise features .feat
4. feature matrix .featmat 
5. results with probability .pairsWprob
6. complexes .cmplx




Adding test line
