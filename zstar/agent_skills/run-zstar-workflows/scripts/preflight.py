#!/usr/bin/env python3
"""Run the packaged ZStar Agent Skill preflight."""

from __future__ import annotations

import argparse

from zstar.agent_skill import DIMENSIONS, LANES, write_preflight_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--lane", choices=LANES, default="bec")
    parser.add_argument("--dim", choices=DIMENSIONS, default="bulk")
    args = parser.parse_args()
    print(write_preflight_json(args.root, lane=args.lane, dimensionality=args.dim))


if __name__ == "__main__":
    main()
