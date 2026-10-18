import argparse
from pathlib import Path


def arguments(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--data-root", type=Path, default=Path("data"))
    parser.add_argument("--output-root", type=Path, default=Path("results"))
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


def result_directory(root, name):
    path = root / name
    path.mkdir(parents=True, exist_ok=True)
    return path
