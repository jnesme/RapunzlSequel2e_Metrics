#!/usr/bin/env python3
"""
PacBio Sequel IIe Run Statistics Parser
Extracts run conditions and statistics from XML files
Handles both CCS/HiFi (ConsensusReadSet) and CLR (SubreadSet) runs
Filters out barcode-specific files to provide run-wide metrics only
"""

import xml.etree.ElementTree as ET
from pathlib import Path
import pandas as pd
from datetime import datetime
import argparse

class SequelRunParser:
    """Parse PacBio Sequel IIe run statistics"""
    
    # XML namespaces
    NAMESPACES = {
        'pbds': 'http://pacificbiosciences.com/PacBioDatasets.xsd',
        'pbbase': 'http://pacificbiosciences.com/PacBioBaseDataModel.xsd',
        'pbmeta': 'http://pacificbiosciences.com/PacBioCollectionMetadata.xsd',
        'ns': 'http://pacificbiosciences.com/PacBioBaseDataModel.xsd'
    }
    
    def __init__(self, base_dir):
        self.base_dir = Path(base_dir)
        
    def find_xml_files(self):
        """Find all relevant XML files in directory structure, excluding barcode-specific files"""
        # Find all XML files recursively
        all_consensus_files = list(self.base_dir.rglob("*.consensusreadset.xml"))
        all_subread_files = list(self.base_dir.rglob("*.subreadset.xml"))
        all_sts_files = list(self.base_dir.rglob("*.sts.xml"))
        
        # Filter out barcode-specific files
        # Barcode files are typically in subdirectories named like "bc2011--bc2011"
        # or have barcode patterns in their filenames
        def is_main_file(filepath):
            """Check if file is a main run file (not barcode-specific)"""
            path_str = str(filepath)
            filename = filepath.name
            
            # Skip if in a barcode subdirectory (e.g., bc2011--bc2011/)
            if '/bc' in path_str and '--bc' in path_str:
                return False
            
            # Skip if filename contains barcode pattern (e.g., .bc2011--bc2011.)
            if '.bc' in filename and '--bc' in filename:
                return False
            
            # Skip "unbarcoded" files as these are also subsets
            if 'unbarcoded' in filename.lower():
                return False
                
            return True
        
        consensus_files = [f for f in all_consensus_files if is_main_file(f)]
        subread_files = [f for f in all_subread_files if is_main_file(f)]
        sts_files = [f for f in all_sts_files if is_main_file(f)]
        
        print(f"Found {len(all_consensus_files)} total consensusreadset files, {len(consensus_files)} main files")
        print(f"Found {len(all_subread_files)} total subreadset files, {len(subread_files)} main files")
        print(f"Found {len(all_sts_files)} total sts files, {len(sts_files)} main files")
        
        return consensus_files, subread_files, sts_files
    
    def parse_consensusreadset(self, xml_path):
        """Extract run conditions from consensusreadset.xml (HiFi/CCS runs)"""
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        data = {
            'xml_file': xml_path.name,
            'run_path': str(xml_path.parent)
        }
        
        # Basic dataset info
        total_length = root.find('.//pbds:TotalLength', self.NAMESPACES)
        num_records = root.find('.//pbds:NumRecords', self.NAMESPACES)
        
        if total_length is not None:
            data['total_length'] = int(total_length.text)
        if num_records is not None:
            data['num_records'] = int(num_records.text)
            
        # Collection metadata
        collection = root.find('.//pbmeta:CollectionMetadata', self.NAMESPACES)
        if collection is not None:
            data['instrument_id'] = collection.get('InstrumentId')
            data['instrument_name'] = collection.get('InstrumentName')
            data['context'] = collection.get('Context')
            data['created_at'] = collection.get('CreatedAt')
            
        # Run details
        run_details = root.find('.//pbmeta:RunDetails', self.NAMESPACES)
        if run_details is not None:
            name = run_details.find('pbmeta:Name', self.NAMESPACES)
            data['run_name'] = name.text if name is not None else None
            
        # Well sample info
        well = root.find('.//pbmeta:WellSample', self.NAMESPACES)
        if well is not None:
            data['sample_name'] = well.get('Name')
            well_name = well.find('pbmeta:WellName', self.NAMESPACES)
            data['well_name'] = well_name.text if well_name is not None else None
            
            insert_size = well.find('pbmeta:InsertSize', self.NAMESPACES)
            data['insert_size'] = int(insert_size.text) if insert_size is not None else None
            
            on_plate_conc = well.find('pbmeta:OnPlateLoadingConcentration', self.NAMESPACES)
            data['loading_concentration'] = float(on_plate_conc.text) if on_plate_conc is not None else None
            
            application = well.find('pbmeta:Application', self.NAMESPACES)
            data['application'] = application.text if application is not None else None
        
        # Automation parameters (Movie length, etc.)
        automation = root.find('.//pbmeta:Automation', self.NAMESPACES)
        if automation is not None:
            for param in automation.findall('.//pbbase:AutomationParameter', self.NAMESPACES):
                param_name = param.get('Name')
                param_value = param.get('SimpleValue')
                
                if param_name == 'MovieLength':
                    data['movie_length_min'] = float(param_value)
                elif param_name == 'InsertSize':
                    data['target_insert_size'] = int(param_value)
                    
        # Kit information
        binding_kit = root.find('.//pbmeta:BindingKit', self.NAMESPACES)
        if binding_kit is not None:
            data['binding_kit'] = binding_kit.get('Name')
            data['binding_kit_part'] = binding_kit.get('PartNumber')
            
        cell_pac = root.find('.//pbmeta:CellPac', self.NAMESPACES)
        if cell_pac is not None:
            data['cell_type'] = cell_pac.get('Name')
            
        return data
    
    def parse_subreadset(self, xml_path):
        """Extract run conditions from subreadset.xml (CLR runs)"""
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        data = {
            'xml_file': xml_path.name,
            'run_path': str(xml_path.parent)
        }
        
        # Basic dataset info
        total_length = root.find('.//pbds:TotalLength', self.NAMESPACES)
        num_records = root.find('.//pbds:NumRecords', self.NAMESPACES)
        
        if total_length is not None:
            data['total_length'] = int(total_length.text)
        if num_records is not None:
            data['num_records'] = int(num_records.text)
            
        # Collection metadata
        collection = root.find('.//pbmeta:CollectionMetadata', self.NAMESPACES)
        if collection is not None:
            data['instrument_id'] = collection.get('InstrumentId')
            data['instrument_name'] = collection.get('InstrumentName')
            data['context'] = collection.get('Context')
            data['created_at'] = collection.get('CreatedAt')
            
        # Run details
        run_details = root.find('.//pbmeta:RunDetails', self.NAMESPACES)
        if run_details is not None:
            name = run_details.find('pbmeta:Name', self.NAMESPACES)
            data['run_name'] = name.text if name is not None else None
            
        # Well sample info
        well = root.find('.//pbmeta:WellSample', self.NAMESPACES)
        if well is not None:
            data['sample_name'] = well.get('Name')
            well_name = well.find('pbmeta:WellName', self.NAMESPACES)
            data['well_name'] = well_name.text if well_name is not None else None
            
            insert_size = well.find('pbmeta:InsertSize', self.NAMESPACES)
            data['insert_size'] = int(insert_size.text) if insert_size is not None else None
            
            on_plate_conc = well.find('pbmeta:OnPlateLoadingConcentration', self.NAMESPACES)
            data['loading_concentration'] = float(on_plate_conc.text) if on_plate_conc is not None else None
        
        # Automation parameters
        automation = root.find('.//pbmeta:Automation', self.NAMESPACES)
        if automation is not None:
            for param in automation.findall('.//pbbase:AutomationParameter', self.NAMESPACES):
                param_name = param.get('Name')
                param_value = param.get('SimpleValue')
                
                if param_name == 'MovieLength':
                    data['movie_length_min'] = float(param_value)
                elif param_name == 'InsertSize':
                    data['target_insert_size'] = int(param_value)
                    
        # Kit information
        binding_kit = root.find('.//pbmeta:BindingKit', self.NAMESPACES)
        if binding_kit is not None:
            data['binding_kit'] = binding_kit.get('Name')
            data['binding_kit_part'] = binding_kit.get('PartNumber')
            
        cell_pac = root.find('.//pbmeta:CellPac', self.NAMESPACES)
        if cell_pac is not None:
            data['cell_type'] = cell_pac.get('Name')
            
        return data
    
    def parse_sts_file(self, xml_path):
        """Extract statistics from .sts.xml file"""
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        data = {
            'sts_file': xml_path.name,
        }
        
        # Helper function to find elements ignoring namespace
        def find_element(parent, tag_name):
            """Find element by tag name, ignoring namespace"""
            for elem in parent:
                # Get tag without namespace
                local_tag = elem.tag.split('}')[-1] if '}' in elem.tag else elem.tag
                if local_tag == tag_name:
                    return elem
            return None
        
        def find_all_elements(parent, tag_name):
            """Find all elements by tag name, ignoring namespace"""
            results = []
            for elem in parent.iter():
                local_tag = elem.tag.split('}')[-1] if '}' in elem.tag else elem.tag
                if local_tag == tag_name:
                    results.append(elem)
            return results
        
        # Basic stats - these are direct children of root
        num_zmws = find_element(root, 'NumSequencingZmws')
        data['num_sequencing_zmws'] = int(num_zmws.text) if num_zmws is not None else None
        
        movie_length = find_element(root, 'MovieLength')
        data['actual_movie_length_min'] = int(movie_length.text) if movie_length is not None else None
        
        sequencing_umy = find_element(root, 'SequencingUmy')
        data['total_bases'] = int(sequencing_umy.text) if sequencing_umy is not None else None
        data['yield_gb'] = round(data['total_bases'] / 1e9, 2) if data['total_bases'] else None
        
        # P0, P1, P2 from LoadingDist
        loading_dist = find_element(root, 'LoadingDist')
        if loading_dist is not None:
            bin_counts = find_all_elements(loading_dist, 'BinCount')
            
            if len(bin_counts) >= 3:
                data['p0_empty'] = int(bin_counts[0].text)
                data['p1_single'] = int(bin_counts[1].text)
                data['p2_multi'] = int(bin_counts[2].text)
                
                if data['num_sequencing_zmws']:
                    data['p0_percent'] = round(100 * data['p0_empty'] / data['num_sequencing_zmws'], 2)
                    data['p1_percent'] = round(100 * data['p1_single'] / data['num_sequencing_zmws'], 2)
                    data['p2_percent'] = round(100 * data['p2_multi'] / data['num_sequencing_zmws'], 2)
        
        # Productivity
        prod_dist = find_element(root, 'ProdDist')
        if prod_dist is not None:
            bin_counts = find_all_elements(prod_dist, 'BinCount')
            if len(bin_counts) >= 2:
                data['productive_zmws'] = int(bin_counts[1].text)  # Index 1 is "Productive"
                if data['num_sequencing_zmws']:
                    data['productivity_percent'] = round(100 * data['productive_zmws'] / data['num_sequencing_zmws'], 2)
        
        # Read length statistics
        read_len_dist = find_element(root, 'ReadLenDist')
        if read_len_dist is not None:
            mean_elem = find_element(read_len_dist, 'SampleMean')
            median_elem = find_element(read_len_dist, 'SampleMed')
            n50_elem = find_element(read_len_dist, 'SampleN50')
            
            data['mean_read_length'] = int(float(mean_elem.text)) if mean_elem is not None else None
            data['median_read_length'] = int(float(median_elem.text)) if median_elem is not None else None
            data['n50_read_length'] = int(float(n50_elem.text)) if n50_elem is not None else None
        
        # SNR values - need to find HqRegionSnrDist elements with Channel attribute
        for channel in ['A', 'C', 'G', 'T']:
            # Find all HqRegionSnrDist elements
            snr_dists = find_all_elements(root, 'HqRegionSnrDist')
            for snr_dist in snr_dists:
                if snr_dist.get('Channel') == channel:
                    mean_elem = find_element(snr_dist, 'SampleMean')
                    if mean_elem is not None:
                        data[f'snr_{channel.lower()}'] = round(float(mean_elem.text), 2)
                    break
        
        return data
    
    def parse_all_runs(self):
        """Parse all runs and combine data"""
        consensus_files, subread_files, sts_files = self.find_xml_files()
        
        all_data = []
        
        # Create mapping of run contexts to files
        run_data = {}
        
        # Parse ConsensusReadSet files (HiFi/CCS runs)
        for consensus_file in consensus_files:
            try:
                data = self.parse_consensusreadset(consensus_file)
                data['run_type'] = 'CCS/HiFi'
                context = data.get('context', consensus_file.stem)
                
                if context not in run_data:
                    run_data[context] = {}
                run_data[context].update(data)
                
            except Exception as e:
                print(f"Error parsing {consensus_file}: {e}")
        
        # Parse SubreadSet files (CLR runs)
        for subread_file in subread_files:
            try:
                data = self.parse_subreadset(subread_file)
                data['run_type'] = 'CLR'
                context = data.get('context', subread_file.stem)
                
                if context not in run_data:
                    run_data[context] = {}
                run_data[context].update(data)
                
            except Exception as e:
                print(f"Error parsing {subread_file}: {e}")
        
        # Match with STS files
        for sts_file in sts_files:
            try:
                data = self.parse_sts_file(sts_file)
                
                # Try to match with consensus data by context
                context = sts_file.stem.split('.')[0]  # e.g., m64241e_251015_113449
                
                if context in run_data:
                    run_data[context].update(data)
                else:
                    run_data[context] = data
                    
            except Exception as e:
                print(f"Error parsing {sts_file}: {e}")
        
        # Convert to list
        all_data = list(run_data.values())
        
        return pd.DataFrame(all_data)

def main():
    parser = argparse.ArgumentParser(
        description='Parse PacBio Sequel IIe run statistics',
        epilog='Example: python parse_sequel_runs.py /path/to/sequel/data -o runs.csv --excel'
    )
    parser.add_argument(
        'base_dir',
        help='Base directory containing Sequel IIe run folders'
    )
    parser.add_argument(
        '-o', '--output',
        default='sequel_runs_summary.csv',
        help='Output CSV file (default: sequel_runs_summary.csv)'
    )
    parser.add_argument(
        '--excel',
        action='store_true',
        help='Also output as Excel file'
    )
    
    args = parser.parse_args()
    
    # Parse runs
    print(f"Scanning {args.base_dir}...")
    print("="*60)
    parser_obj = SequelRunParser(args.base_dir)
    df = parser_obj.parse_all_runs()
    
    if len(df) == 0:
        print("\nNo runs found!")
        return
    
    # Sort by date
    if 'created_at' in df.columns:
        df['created_at'] = pd.to_datetime(df['created_at'], errors='coerce')
        df = df.sort_values('created_at', ascending=False)
    
    # Save output
    df.to_csv(args.output, index=False)
    print(f"\n✓ Saved {len(df)} runs to {args.output}")
    
    if args.excel:
        excel_output = args.output.replace('.csv', '.xlsx')
        df.to_excel(excel_output, index=False, engine='openpyxl')
        print(f"✓ Saved Excel output to {excel_output}")
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total runs parsed: {len(df)}")
    
    if 'run_type' in df.columns:
        print("\nRun types:")
        print(df['run_type'].value_counts().to_string())
    
    if 'yield_gb' in df.columns:
        total_yield = df['yield_gb'].sum()
        avg_yield = df['yield_gb'].mean()
        print(f"\nTotal yield: {total_yield:.1f} Gb")
        print(f"Average yield per run: {avg_yield:.1f} Gb")
    
    if 'p1_percent' in df.columns:
        avg_p1 = df['p1_percent'].mean()
        print(f"\nAverage P1 (single loading): {avg_p1:.1f}%")
    
    if 'productivity_percent' in df.columns:
        avg_prod = df['productivity_percent'].mean()
        print(f"Average productivity: {avg_prod:.1f}%")
    
    print("\n" + "="*60)
    print(f"Columns extracted ({len(df.columns)}):")
    print("="*60)
    for col in sorted(df.columns):
        non_null = df[col].notna().sum()
        print(f"  • {col:30s} ({non_null}/{len(df)} populated)")

if __name__ == '__main__':
    main()
