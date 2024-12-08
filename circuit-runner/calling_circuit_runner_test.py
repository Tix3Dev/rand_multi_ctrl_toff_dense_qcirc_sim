import subprocess

def run_command():
    command = ["./target/release/main", "python_generated_dot.txt"]
    result = subprocess.run(command, capture_output=True, text=True)
    print(result.stdout)

if __name__ == "__main__":
    run_command()