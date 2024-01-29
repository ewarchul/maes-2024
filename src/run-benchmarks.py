import subprocess
import os
import shlex
import argparse
from pathlib import Path

SOURCE_DIR = Path(os.path.dirname(os.path.abspath(__file__)))

def spawn_benchmark_runner(benchmark_spec: Path): 
    args = shlex.split(f'R --no-save -e "library(cecb); cecb::run_benchmark(\\"{benchmark_spec}\\")"')
    subprocess.Popen(args) 

def distpach_machine(machine_name: str) -> Path:
    if machine_name == "hal":
        return SOURCE_DIR / "yml" / "hal"
    else:
        return SOURCE_DIR / "yml" / "borg"



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="cec_runner",
        description="starts CEC benchmark using yaml specification"
    )
    parser.add_argument('-m', '--machine', choices=['hal', 'borg'], required=True)
    args = parser.parse_args()

    for spec in distpach_machine(args.machine).iterdir():
        spawn_benchmark_runner(spec)



