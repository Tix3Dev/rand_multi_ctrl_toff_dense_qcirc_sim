# How to benchmark (for now)

1. Open `my_had_pyzx_to_quizx_dot_repr.ipynb` and `poc.ipynb` (CNOT)

2. In both, initialize the same circuit, e.g.:
```python
my_seed = random.randint(1, 100000000)
my_seed = 77900050 # If you like what you see, copy the generated seed an paste it here, so next time you run it, it will be the same
print("Randomly generated seed:", my_seed)

# randqcirc = RandQCirc(qbit_cnt=7, cnot_cnt=20, tof_cnt=5, allow_cnot_in_targ_range=False, tof_targ_qbit_range=range(5,6+1), allow_ctrl_in_targ_range=False, seed=my_seed)
randqcirc = RandQCirc(qbit_cnt=100, cnot_cnt=200, tof_cnt=50, allow_cnot_in_targ_range=False, tof_targ_qbit_range=range(90,99+1), allow_ctrl_in_targ_range=False, seed=my_seed)
randqcirc.generate_string_circuit()
randqcirc.print_string_circ_repr()
randqcirc.generate_diagram()

# randqcirc.g.apply_state("+"*5 + "0"*2)
# randqcirc.g.apply_effect("+"*5 + "0"*2)
randqcirc.g.apply_state("+"*90 + "0"*10)
randqcirc.g.apply_effect("+"*90 + "0"*10)
```

3. Run `poc.ipynb`, at the bottom, you will see how long it took and how many terms got produced

4. Run `my_had_pyzx_to_quizx_dot_repr.ipynb`. This will generate `python_generated_dot.txt`.

5. In `circuit-runner/` (Ubuntu WSL), run `cargo build --releae` if not already done. This step only has to be done once

6. In `circuit-runner/` (Ubuntu WSL), run `./target/release/main python_generated_dot.txt`. At the bottom, you will see how long it took and how many terms got produced

7. Compare 3. and 7.