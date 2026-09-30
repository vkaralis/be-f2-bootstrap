"""Command-line interface for headerless numeric CSV inputs."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from .core import bootstrap_f2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--times", type=float, nargs="+", required=True, help="actual shared sampling times in minutes")
    parser.add_argument("--framework", choices=["ema", "ich-m13b"], required=True)
    parser.add_argument("--num-bootstraps", type=int, default=30000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    try:
        test = np.loadtxt(args.test, delimiter=",", ndmin=2)
        reference = np.loadtxt(args.reference, delimiter=",", ndmin=2)
        result = bootstrap_f2(test, reference, args.times, framework=args.framework,
                              num_bootstraps=args.num_bootstraps, seed=args.seed)
        args.output.mkdir(parents=True, exist_ok=True)
        np.savetxt(args.output / "test.csv", test, delimiter=",", fmt="%.17g")
        np.savetxt(args.output / "reference.csv", reference, delimiter=",", fmt="%.17g")
        (args.output / "times.json").write_text(json.dumps(args.times) + "\n", encoding="utf-8")
        (args.output / "analysis.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        with (args.output / "bootstrap_f2.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["f2"])
            writer.writerows([value] for value in result["bootstrap_f2"])
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps({"ci90": result["ci90"], "bootstrap_criterion_met": result["bootstrap_criterion_met"]}))


if __name__ == "__main__":
    main()
