"""Expand SVG ``use`` elements for reliable PowerPoint shape conversion."""

from __future__ import annotations

import argparse
import copy
from pathlib import Path
import xml.etree.ElementTree as ET


SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
USE_TAG = f"{{{SVG_NS}}}use"


def expand_uses(source: Path, destination: Path) -> int:
    ET.register_namespace("", SVG_NS)
    ET.register_namespace("xlink", XLINK_NS)
    tree = ET.parse(source)
    root = tree.getroot()
    references = {
        element_id: element
        for element in root.iter()
        if (element_id := element.get("id"))
    }
    expanded = 0

    while True:
        parent_map = {child: parent for parent in root.iter() for child in parent}
        uses = list(root.iter(USE_TAG))
        if not uses:
            break
        for use in uses:
            href = use.get(f"{{{XLINK_NS}}}href") or use.get("href")
            if not href or not href.startswith("#") or href[1:] not in references:
                raise ValueError(f"Unresolved SVG use reference in {source}: {href!r}")
            target = copy.deepcopy(references[href[1:]])
            target.attrib.pop("id", None)

            transforms = []
            if transform := use.get("transform"):
                transforms.append(transform)
            x = use.get("x", "0")
            y = use.get("y", "0")
            if x != "0" or y != "0":
                transforms.append(f"translate({x} {y})")
            if transform := target.get("transform"):
                transforms.append(transform)
            if transforms:
                target.set("transform", " ".join(transforms))
            else:
                target.attrib.pop("transform", None)

            ignored = {"href", f"{{{XLINK_NS}}}href", "x", "y", "transform"}
            for name, value in use.attrib.items():
                if name not in ignored:
                    target.set(name, value)

            parent = parent_map[use]
            parent.insert(list(parent).index(use), target)
            parent.remove(use)
            expanded += 1

    destination.parent.mkdir(parents=True, exist_ok=True)
    tree.write(destination, encoding="utf-8", xml_declaration=True)
    return expanded


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    expanded = expand_uses(args.source.resolve(), args.destination.resolve())
    print(f"Expanded {expanded} SVG use elements: {args.destination}")


if __name__ == "__main__":
    main()
