#!/usr/bin/env python3
"""Check for different XML structures across .sts.xml files"""

import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict

base_dir = Path('/projects/codon_0000/data/RapunzlSequel2e/')
sts_files = list(base_dir.rglob("*.sts.xml"))

# Filter out barcode files
def is_main_file(filepath):
    path_str = str(filepath)
    filename = filepath.name
    if '/bc' in path_str and '--bc' in path_str:
        return False
    if '.bc' in filename and '--bc' in filename:
        return False
    if 'unbarcoded' in filename.lower():
        return False
    return True

sts_files = [f for f in sts_files if is_main_file(f)]

print(f"Analyzing {len(sts_files)} main .sts.xml files...")
print()

# Group files by XML structure
structures = defaultdict(list)

for sts_file in sts_files[:20]:  # Check first 20 files
    try:
        tree = ET.parse(sts_file)
        root = tree.getroot()
        
        # Get version if available
        version = root.get('Version', 'unknown')
        
        # Get root tag (without namespace)
        root_tag = root.tag.split('}')[-1] if '}' in root.tag else root.tag
        
        # Check for key elements
        has_num_zmws = False
        has_movie_length = False
        has_sequencing_umy = False
        
        for child in root:
            tag = child.tag.split('}')[-1] if '}' in child.tag else child.tag
            if tag == 'NumSequencingZmws':
                has_num_zmws = True
            elif tag == 'MovieLength':
                has_movie_length = True
            elif tag == 'SequencingUmy':
                has_sequencing_umy = True
        
        structure_key = f"{root_tag}_v{version}_zmws{has_num_zmws}_movie{has_movie_length}_umy{has_sequencing_umy}"
        structures[structure_key].append(sts_file)
        
    except Exception as e:
        print(f"Error parsing {sts_file.name}: {e}")

print(f"Found {len(structures)} different XML structures:")
print()

for structure, files in structures.items():
    print(f"Structure: {structure}")
    print(f"  Files: {len(files)}")
    print(f"  Example: {files[0].name}")
    
    # Parse one example to show details
    try:
        tree = ET.parse(files[0])
        root = tree.getroot()
        
        # Show first few direct children
        print(f"  Root children (first 10):")
        for i, child in enumerate(root):
            if i >= 10:
                break
            tag = child.tag.split('}')[-1] if '}' in child.tag else child.tag
            text = child.text[:50] if child.text and len(child.text) > 50 else child.text
            print(f"    {tag}: {text}")
    except:
        pass
    
    print()
