#!/usr/bin/env python3
"""Debug script to check .sts.xml file parsing"""

import sys
from pathlib import Path

# Remove old module if cached
if 'parse_sequel_runs' in sys.modules:
    del sys.modules['parse_sequel_runs']

sys.path.insert(0, '/mnt/user-data/outputs')

from parse_sequel_runs import SequelRunParser

# Initialize parser
parser = SequelRunParser('/projects/codon_0000/data/RapunzlSequel2e/')

# Find files
consensus_files, subread_files, sts_files = parser.find_xml_files()

print(f"Found {len(sts_files)} .sts.xml files")
print()

if sts_files:
    print("Testing first .sts.xml file...")
    sts_file = sts_files[0]
    print(f"File: {sts_file.name}")
    print()
    
    # Try parsing
    if sts_file.exists():
        try:
            data = parser.parse_sts_file(sts_file)
            print("✓ Parsed successfully!")
            print()
            print("Extracted data:")
            for key, value in sorted(data.items()):
                if value is not None:
                    print(f"  ✓ {key:30s} = {value}")
                else:
                    print(f"  ✗ {key:30s} = None")
        except Exception as e:
            print(f"✗ ERROR: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("ERROR: File doesn't exist!")
