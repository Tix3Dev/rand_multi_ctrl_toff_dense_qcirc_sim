# How to run cherry picked examples using quizx

1. Generate `python_generated_dot.txt` using `cherry_pick_analysis.ipynb`.

2. In `circuit-runner/` (Ubuntu WSL), run `cargo build --releae` if not already done. This step only has to be done once

3. In `circuit-runner/` (Ubuntu WSL), run `./target/release/main python_generated_dot.txt`. At the bottom, you will see how long it took and how many terms got produced
