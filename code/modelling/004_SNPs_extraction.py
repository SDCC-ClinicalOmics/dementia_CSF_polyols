import pandas as pd
import numpy as np
import os
from glob import glob
from tqdm import tqdm

# ----------------------------------
# Configuration
# ----------------------------------
ILLUMINA_DIR = "ADNI2"
OUTPUT_DIR = "."
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Primary candidate SNPs (what we ideally want)
PRIMARY_SNPS = {
    'rs2229540': {'effect_allele': 'G', 'chr': 1, 'pos': 46032311, 'gene': 'AKR1A1', 'metabolite': 'Erythritol'},
    'rs7239738': {'effect_allele': 'A', 'chr': 18, 'pos': 680520, 'gene': 'ENOSF1', 'metabolite': 'Ribonic acid'},
    'rs3786349': {'effect_allele': 'A', 'chr': 18, 'pos': 712568, 'gene': 'ENOSF1', 'metabolite': 'Ribonic acid'},
    'rs55901542': {'effect_allele': 'C', 'chr': 15, 'pos': 45069019, 'gene': 'SORD', 'metabolite': 'Arabitol/Ribitol'},
    'rs2790': {'effect_allele': 'A', 'chr': 18, 'pos': 673086, 'gene': 'ENOSF1', 'metabolite': 'Ribonic acid'},
    'rs2298581': {'effect_allele': 'G', 'chr': 18, 'pos': 677931, 'gene': 'ENOSF1', 'metabolite': 'Ribonic acid'},
    'rs10502289': {'effect_allele': 'A', 'chr': 18, 'pos': 676789, 'gene': 'ENOSF1', 'metabolite': 'Ribonic acid'},
    'rs2847153': {'effect_allele': 'A', 'chr': 18, 'pos': 661647, 'gene': 'TYMS', 'metabolite': 'Arabonate'},
    'rs377697486': {'effect_allele': 'C', 'chr': 15, 'pos': 45046133, 'gene': 'SORD', 'metabolite': 'Ribitol'},
    'rs148563337': {'effect_allele': 'T', 'chr': 1, 'pos': 45559808, 'gene': 'AKR1A1', 'metabolite': 'Erythritol'},
}

# ----------------------------------
# Functions
# ----------------------------------
def genotype_to_additive(allele1, allele2, effect_allele):
    """Convert genotype to additive coding (0, 1, 2)"""
    if pd.isna(allele1) or pd.isna(allele2):
        return np.nan
    if str(allele1) in ['0', '-', ''] or str(allele2) in ['0', '-', '']:
        return np.nan
    count = int(allele1 == effect_allele) + int(allele2 == effect_allele)
    return count

def find_available_and_proxies(sample_file, primary_snps, proxy_window=50000):
    """
    Check which primary SNPs exist, and find proxies for missing ones
    """
    print(f"Reading sample file: {os.path.basename(sample_file)}")
    df = pd.read_csv(sample_file)
    print(f"Total SNPs in file: {len(df)}")
    
    available = {}
    proxies = {}
    
    # Get list of all SNP names in the file for quick lookup
    available_snp_names = set(df['SNP Name'].values)
    
    for rsid, info in primary_snps.items():
        # Check if primary SNP exists
        if rsid in available_snp_names:
            available[rsid] = info
            print(f"  ✓ Found {rsid}")
        else:
            print(f"  ✗ Missing {rsid} - searching for proxies...")
            # Find proxy SNPs in the region
            chr_num = info['chr']
            pos = info['pos']
            
            region = df[(df['Chr'] == chr_num) & 
                       (df['Position'] >= pos - proxy_window) & 
                       (df['Position'] <= pos + proxy_window)].copy()
            
            print(f"    Region chr{chr_num}:{pos-proxy_window:,}-{pos+proxy_window:,}: {len(region)} SNPs")
            
            if len(region) > 0:
                region['distance'] = abs(region['Position'] - pos)
                region = region.sort_values('distance')
                
                # Get top 5 closest SNPs as potential proxies
                proxy_candidates = []
                for _, row in region.head(5).iterrows():
                    proxy_candidates.append({
                        'rsid': row['SNP Name'],
                        'position': int(row['Position']),
                        'distance': int(row['distance']),
                        'alleles': row['SNP'],
                        'original_snp': rsid,
                        'gene': info['gene'],
                        'metabolite': info['metabolite']
                    })
                
                proxies[rsid] = proxy_candidates
                print(f"    Found {len(proxy_candidates)} proxy candidates")
            else:
                proxies[rsid] = []
                print(f"    ⚠ No SNPs found in region")
    
    return available, proxies

def extract_snps_from_file(filepath, snps_to_extract):
    """Extract SNPs from a single Illumina CSV file"""
    try:
        df = pd.read_csv(filepath)
        sample_id = df['Sample ID'].iloc[0]
        
        df_candidates = df[df['SNP Name'].isin(snps_to_extract.keys())]
        
        snp_scores = {}
        for _, row in df_candidates.iterrows():
            snp_name = row['SNP Name']
            effect_allele = snps_to_extract[snp_name]
            
            allele1 = row['Allele1 - Forward']
            allele2 = row['Allele2 - Forward']
            
            score = genotype_to_additive(allele1, allele2, effect_allele)
            snp_scores[snp_name] = score
        
        return sample_id, snp_scores
    
    except Exception as e:
        return None, None

def process_all_files(illumina_dir, snps_to_extract):
    """Process all Illumina CSV files"""
    file_pattern = os.path.join(illumina_dir, "*.csv")
    files = glob(file_pattern)
    
    print(f"\nFound {len(files)} CSV files to process")
    
    if len(files) == 0:
        return pd.DataFrame()
    
    all_samples = {}
    failed_files = 0
    
    for filepath in tqdm(files, desc="Processing samples"):
        sample_id, snp_scores = extract_snps_from_file(filepath, snps_to_extract)
        
        if sample_id is not None:
            all_samples[sample_id] = snp_scores
        else:
            failed_files += 1
    
    if failed_files > 0:
        print(f"Warning: {failed_files} files failed")
    
    df = pd.DataFrame.from_dict(all_samples, orient='index')
    return df

# ----------------------------------
# Main Processing
# ----------------------------------
if __name__ == "__main__":
    
    print("="*80)
    print("INTEGRATED SNP EXTRACTION + PROXY FINDING")
    print("="*80)
    
    # Step 1: Find available SNPs and proxies
    print("\n" + "="*80)
    print("STEP 1: SCANNING FOR AVAILABLE SNPs AND PROXIES")
    print("="*80 + "\n")
    
    sample_files = glob(os.path.join(ILLUMINA_DIR, "*.csv"))
    
    if len(sample_files) == 0:
        print(f"ERROR: No CSV files found in {ILLUMINA_DIR}")
        exit(1)
    
    test_file = sample_files[0]
    available, proxies = find_available_and_proxies(test_file, PRIMARY_SNPS, proxy_window=50000)
    
    # Report findings
    print(f"\n{'='*80}")
    print(f"SUMMARY: {len(available)}/{len(PRIMARY_SNPS)} SNPs directly available")
    print(f"{'='*80}")
    
    print(f"\n✓ DIRECTLY AVAILABLE ({len(available)} SNPs):")
    for rsid, info in available.items():
        print(f"  {rsid:15s} - {info['gene']:10s} -> {info['metabolite']}")
    
    print(f"\n✗ MISSING ({len(proxies)} SNPs) - PROXIES:")
    for missing_rsid, candidates in proxies.items():
        info = PRIMARY_SNPS[missing_rsid]
        print(f"\n  {missing_rsid} ({info['gene']} -> {info['metabolite']}, chr{info['chr']}:{info['pos']:,}):")
        if candidates:
            print(f"    Proxy candidates:")
            for i, proxy in enumerate(candidates[:3], 1):
                print(f"      {i}. {proxy['rsid']:15s} distance={proxy['distance']:>7,} bp  alleles={proxy['alleles']}")
        else:
            print(f"    ⚠ NO proxies found within ±50kb")
    
    # Step 2: Build final SNP set
    print(f"\n{'='*80}")
    print("STEP 2: BUILDING FINAL SNP SET")
    print(f"{'='*80}\n")
    
    final_snps = {}
    proxy_mapping = {}
    
    # Add directly available SNPs
    for rsid, info in available.items():
        final_snps[rsid] = info['effect_allele']
        print(f"  ✓ {rsid:15s} (direct, effect allele = {info['effect_allele']})")
    
    # Add best proxy for each missing SNP
    for missing_rsid, candidates in proxies.items():
        if candidates and len(candidates) > 0:
            best_proxy = candidates[0]
            proxy_rsid = best_proxy['rsid']
            
            # Parse alleles from the [A/G] format
            alleles_str = best_proxy['alleles'].strip('[]').split('/')
            
            # For now, use the first allele (you'll need to verify direction)
            proxy_effect_allele = alleles_str[0]
            
            final_snps[proxy_rsid] = proxy_effect_allele
            proxy_mapping[proxy_rsid] = {
                'proxy_for': missing_rsid,
                'distance': best_proxy['distance'],
                'gene': best_proxy['gene'],
                'metabolite': best_proxy['metabolite'],
                'position': best_proxy['position'],
                'alleles': best_proxy['alleles'],
                'selected_effect_allele': proxy_effect_allele
            }
            print(f"  ⊕ {proxy_rsid:15s} (proxy for {missing_rsid}, dist={best_proxy['distance']:,} bp, alleles={best_proxy['alleles']}, using {proxy_effect_allele})")
        else:
            print(f"  ✗ {missing_rsid:15s} (no proxies available - EXCLUDED)")
    
    print(f"\nFinal SNP set: {len(final_snps)} SNPs ({len(available)} direct + {len(proxy_mapping)} proxies)")
    
    # Save proxy mapping
    if proxy_mapping:
        proxy_df = pd.DataFrame.from_dict(proxy_mapping, orient='index')
        proxy_file = os.path.join(OUTPUT_DIR, "proxy_snp_mapping.csv")
        proxy_df.to_csv(proxy_file)
        print(f"\n⚠ IMPORTANT: Proxy effect alleles auto-selected - VERIFY in: {proxy_file}")
    
    # Step 3: Extract genotypes
    print(f"\n{'='*80}")
    print("STEP 3: EXTRACTING GENOTYPES FROM ALL SAMPLES")
    print(f"{'='*80}")
    
    snp_data = process_all_files(ILLUMINA_DIR, final_snps)
    
    if snp_data.empty:
        print("ERROR: No data extracted")
        exit(1)
    
    # QC
    print(f"\n{'='*80}")
    print(f"EXTRACTION COMPLETE: {len(snp_data)} samples × {len(snp_data.columns)} SNPs")
    print(f"{'='*80}")
    
    print("\nMissingness by SNP:")
    for snp in snp_data.columns:
        missing_count = snp_data[snp].isnull().sum()
        missing_pct = (missing_count / len(snp_data)) * 100
        status = "PROXY" if snp in proxy_mapping else "DIRECT"
        print(f"  {snp:15s} ({status:6s}): {missing_count:4d} ({missing_pct:5.1f}%)")
    
    print("\nAllele frequencies:")
    for snp in snp_data.columns:
        valid = snp_data[snp].dropna()
        if len(valid) > 0:
            freq = valid.mean() / 2
            status = "PROXY" if snp in proxy_mapping else "DIRECT"
            geno_dist = valid.value_counts().sort_index()
            dist_str = " ".join([f"{int(g)}:{int(c)}" for g, c in geno_dist.items()])
            print(f"  {snp:15s} ({status:6s}): EAF={freq:.3f}  [{dist_str}]")
    
    # Save outputs
    print(f"\n{'='*80}")
    print("SAVING OUTPUTS")
    print(f"{'='*80}")
    
    snp_data.to_csv(os.path.join(OUTPUT_DIR, "candidate_snps_genotypes.csv"))
    print("  ✓ candidate_snps_genotypes.csv")
    
    snp_data_reset = snp_data.reset_index().rename(columns={'index': 'Sample_ID'})
    snp_data_reset.to_csv(os.path.join(OUTPUT_DIR, "candidate_snps_for_merge.csv"), index=False)
    print("  ✓ candidate_snps_for_merge.csv")
    
    if proxy_mapping:
        print("  ✓ proxy_snp_mapping.csv")
    
    # Metadata
    metadata = []
    for rsid in snp_data.columns:
        if rsid in available:
            metadata.append({
                'SNP': rsid,
                'Type': 'Direct',
                'Gene': available[rsid]['gene'],
                'Chr': available[rsid]['chr'],
                'Position': available[rsid]['pos'],
                'Metabolite': available[rsid]['metabolite'],
                'Effect_Allele': final_snps[rsid],
                'Proxy_For': '',
                'Proxy_Distance': ''
            })
        elif rsid in proxy_mapping:
            pm = proxy_mapping[rsid]
            metadata.append({
                'SNP': rsid,
                'Type': 'Proxy',
                'Gene': pm['gene'],
                'Chr': PRIMARY_SNPS[pm['proxy_for']]['chr'],
                'Position': pm['position'],
                'Metabolite': pm['metabolite'],
                'Effect_Allele': pm['selected_effect_allele'],
                'Proxy_For': pm['proxy_for'],
                'Proxy_Distance': pm['distance']
            })
    
    metadata_df = pd.DataFrame(metadata)
    metadata_df.to_csv(os.path.join(OUTPUT_DIR, "snp_metadata.csv"), index=False)
    print("  ✓ snp_metadata.csv")
    
    print(f"\n{'='*80}")
    print("✅ COMPLETE!")
    print(f"{'='*80}")
    print(f"\nYou now have {len(final_snps)} SNPs:")
    print(f"  • {len(available)} direct matches")
    print(f"  • {len(proxy_mapping)} proxy SNPs")
    print("\nNext: Review proxy_snp_mapping.csv and verify effect alleles")