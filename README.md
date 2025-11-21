# Acknowledgements
Professor Brian Hie wrote and applied the ESM scripts. 

# RubiscoMLDE
This repository contains scripts relating to the ML-guided directed evolution of plant Rubisco.

# ESM Code
Code is in the [`bin/`](bin/) directory.
ESM-1b/v refers to the sequence-alone language models, while ESM-IF1 refers to the structure-aware model.

# Results
Recommended mutations for for the sequence-alone model are in [`target/results_esm1bv.txt`](target/results_esm1bv.txt).
Code for producing these results is in [`bin/esm1bv.py`](bin/esm1bv.py).

Recommended mutations for for the structure-aware model are in [`target/results_esmif1_complex_better_1RLC.txt`](target/results_esmif1_complex_better_1RLC.txt) and [`target/results_esmif1_complex_better_4RUB.txt`](target/results_esmif1_complex_better_4RUB.txt).
Code for producing these results is in [`bin/esmif1_score_complex.py`](bin/esmif1_score_complex.py).

# Figure generation
Code and inputs for heatmap plots are in the [`figures/`](figures/) directory.

# Sequence alignments
Code and output alignments is in the [`alignments/`](alignments/) directory.
