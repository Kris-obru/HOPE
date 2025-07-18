import requests
import math
# -------------------------------------------------------------------------------------------------
# FUNCTIONS MODULE FOR GENOMIC MAPPING PIPELINE
# -------------------------------------------------------------------------------------------------
# This file contains all utility and pipeline functions used in the main mapping and analysis pipeline.
# Each section is clearly marked by phase, stage, and step, mirroring the structure in main.py.
# -------------------------------------------------------------------------------------------------

__all__ = [
    'calculate_rarity_modifier',
    'get_snp_types',
]

# -------------------------------------------------------------------------------------------------
# MAP CREATION PHASE
# -------------------------------------------------------------------------------------------------

# ================================================================================================
# Stage 2: Initial Mapping Phase
# ================================================================================================
# ------------------------ Genomic System --------------------------------------------------------
#   - Step 1: SNP annotation and weight calculation
#   - Step 2: Downstream interaction gene retrieval (see string_api.py)
#   - Step 3: (future) Tissue relevance, pathway, and node/edge construction
# -----------------------------------------------------------------------------------------------


def calculate_rarity_modifier(maf, steepness=30, midpoint=0.05, max_boost=1.5):
    """
    Calculates a smooth, capped rarity modifier that boosts rare SNPs.
    Used in: Stage 2 - Genomic system, for calculating rarity modifier for SNP node weights.
    Called from: main.py, when building SNP node attributes.
    
    Parameters:
    - maf: Minor Allele Frequency (0 to 0.5)
    - steepness: how quickly the curve rises as maf decreases
    - midpoint: maf value where modifier = midpoint of range (default 0.05)
    - max_boost: maximum boost allowed (default 1.5)

    Returns:
    - rarity_modifier: multiplier between 1.0 and max_boost
    """
    if maf >= 0.5 or maf < 0:
        return 1.0  # invalid input or common variant
    logistic = 1 / (1 + math.exp(steepness * (maf - midpoint)))
    return 1.0 + (max_boost - 1.0) * logistic

# -----------------------------------------------------------------------------------------------
#   SNP Annotation and Weight Calculation
# -----------------------------------------------------------------------------------------------

def get_snp_types(snp_list, base_weight=1.0):
    """
    Processes multiple SNPs using Ensembl VEP and calculates snp_weights.
    Used in: Stage 2 - Genomic system, for retrieving SNP annotation and calculating SNP weights.
    Called from: main.py, to process the user's SNP list and get gene/consequence/weight info.

    Args:
        snp_list (dict): Dict of { 'rsID': zygosity_weight }
        base_weight (float): Optional scaling weight

    Returns:
        List of dicts, each with:
            'snp_id', 'gene', 'consequence', 'snp_type_weight',
            'zygosity_weight', 'base_weight', 'final_snp_weight'
    """
    SNP_WEIGHTS = {
        "regulatory_region_variant": 1.0,
        "TF_binding_site_variant": 1.0,
        "5_prime_UTR_variant": 0.9,
        "missense_variant": 0.9,
        "splice_region_variant": 0.8,
        "synonymous_variant": 0.4,
        "intron_variant": 0.2,
        "intergenic_variant": 0.1
    }
    # Prepare POST request payload
    vep_url = "https://rest.ensembl.org/vep/human/id"
    headers = { "Content-Type": "application/json" }
    payload = { "ids": list(snp_list.keys()) }

    response = requests.post(vep_url, headers=headers, json=payload)

    if response.status_code != 200:
        raise Exception(f"VEP batch query failed: {response.status_code}")

    vep_data = response.json()
    results = []

    for record in vep_data:
        snp_id = record.get('id', 'UNKNOWN')
        consequence = record.get('most_severe_consequence', 'UNKNOWN')
        transcript_data = record.get('transcript_consequences', [])
        gene_symbol = transcript_data[0].get('gene_symbol', 'UNKNOWN') if transcript_data else 'UNKNOWN'

        snp_type_weight = SNP_WEIGHTS.get(consequence, 0)
        zygosity_weight = snp_list.get(snp_id, 1.0)

        final_snp_weight = base_weight * zygosity_weight * snp_type_weight

        results.append({
            'snp_id': snp_id,
            'gene': gene_symbol,
            'consequence': consequence,
            'snp_type_weight': snp_type_weight,
            'zygosity_weight': zygosity_weight,
            'base_weight': base_weight,
            'final_snp_weight': final_snp_weight
        })

    return results

# -----------------------------------------------------------------------------------------------
# (Future) Additional stages and steps can be added here, following the same commenting structure
# -----------------------------------------------------------------------------------------------
