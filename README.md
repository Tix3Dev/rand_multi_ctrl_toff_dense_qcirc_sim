# rand_multi_ctrl_toff_dense_qcirc_sim

This repository is part of my [Maturitätsarbeit](https://github.com/Tix3Dev/Maturitaetsarbeit).

It contains a [tool](https://github.com/Tix3Dev/rand_multi_ctrl_toff_dense_qcirc_sim/blob/main/circuit-runner/runner.py) to execute a parallelized benchmark and the
corresponding [Jupyter notebook](https://github.com/Tix3Dev/rand_multi_ctrl_toff_dense_qcirc_sim/blob/main/circuit-runner/benchmark_eval.ipynb) to evaluate/plot the results.
Additionally, it contains a [Jupyter notebook](https://github.com/Tix3Dev/rand_multi_ctrl_toff_dense_qcirc_sim/blob/main/circuit-runner/cherry_pick_analysis.ipynb) to test
cherry-picked examples, which can be compared to quizx according to the [instructions](https://github.com/Tix3Dev/rand_multi_ctrl_toff_dense_qcirc_sim/blob/main/circuit-runner/instructions.md).

The aim of the benchmark is to test a slightly refined version of our [multiloop Feynman diagram simulation algorithm](https://github.com/Tix3Dev/feynman_loop_diagram_qsim), which is now applied on randomly generated multi-control
Toffoli gate dense quantum circuits.
