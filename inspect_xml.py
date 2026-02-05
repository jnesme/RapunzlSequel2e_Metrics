#!/usr/bin/env python3
"""Inspect XML structure of .sts.xml file"""

import xml.etree.ElementTree as ET
from pathlib import Path
from collections import Counter

sts_file = Path('/projects/codon_0000/data/RapunzlSequel2e/r64241e_20240314_114036/3_C01/m64241e_240316_184724.sts.xml')

if not sts_file.exists():
    print(f"File not found: {sts_file}")
    exit(1)

print(f"Analyzing: {sts_file.name}")
print(f"Size: {sts_file.stat().st_size} bytes")
print()

tree = ET.parse(sts_file)
root = tree.getroot()

print(f"Root tag: {root.tag}")
print(f"Root attributes: {root.attrib}")
print()

# Collect all tags
all_tags = []
def collect_tags(element, path=""):
    tag = element.tag.split('}')[-1] if '}' in element.tag else element.tag
    current_path = f"{path}/{tag}" if path else tag
    all_tags.append((current_path, tag))
    for child in element:
        collect_tags(child, current_path)

collect_tags(root)

# Count unique tags
tag_counter = Counter([t[1] for t in all_tags])

print("All unique tags in file:")
for tag, count in sorted(tag_counter.items()):
    print(f"  {tag}: {count} occurrences")

print()
print("Looking for expected tags:")
expected = ['NumSequencingZmws', 'MovieLength', 'SequencingUmy', 'LoadingDist', 
            'ProdDist', 'ReadLenDist', 'HqRegionSnrDist']
for exp in expected:
    found = [p for p, t in all_tags if exp in t]
    if found:
        print(f"  ✓ {exp}: Found at {found[0]}")
    else:
        print(f"  ✗ {exp}: NOT FOUND")

print()
print("First 10 paths with interesting names:")
interesting = [p for p, t in all_tags if any(x in p for x in ['Zmw', 'Movie', 'Length', 'Base', 'Loading', 'Prod'])]
for path in interesting[:10]:
    print(f"  {path}")
