import requests
import math
from string_api import StringAPI
from functions import *

# SNP list contains rsIDs and zygosity weights. 0.7 if heterozygous (for functional bias), 1 if homozygous
snp_list = {'rs4680': 0.7, 'rs6323': 1.0, 'rs1801133': 1.0}

# ------------------ MAP CREATION PHASE ------------------

# Stage 1: Retrieve all requisite user information
# - Genetic Data - WGS or clinical grade SNP chip only for accuracy
# - Nutrient & Lab Data
# - Symptom & Health Profile
#   - Detailed symptom log: Mood, energy, cognition, digestion, sleep
#   - Condition history: Diagnosed or suspected conditions
#   - Medication/supplement use: Doses, duration, effects
#   - Environmental exposures: Mold, allergens, stressors
#
# - Lifestyle and Demographic Data
#   - Diet: Macronutrient ratios, allergies, intolerances
#   - Sleep: Quality, duration, timing
#   - Exercise: Type, intensity, frequency
#   - Demographics: Age, sex, ancestry/ethnicity
#
# - Time-Series or Longitudinal Data
#   - Supplement/symptom response timelines
#   - Before/after intervention data
#   - Daily or weekly logs for AI training



######################################## Stage 2: Initial mapping phase #####################################################
#
# ------------------------ Genomic system: ---------------------------------------------------------------------------------
#   - Basic SNP node map (user snps) (label as main snps)
#   - Downstream interaction nodes (any genes that interact with the main snp genes that are above .4 or .2 confidence threshold when multiplied with the connecting gene for hypothesis generation) (label as downstream interaction)
#   - Important genes that interact with the lab data, diagnosed conditions, medications, environmental exposures (label as interaction gene)
#
# --- Step 1:
#    
#   - Use VEP to get SNP type for all main SNPs (There are different types of SNPs, and some have more weight on the system)
#   - Calculate the snp weight.
#   
#   - (needs more comprehensive formula) snp_weight = base_weight × zygosity_weight × snp_type_weight x rarity_weight x functional_weight

snp_node_weight_example = {
    "gene": "MAOA",                # Gene symbol (e.g., "MAOA")
    "snp_weight": 0.4,              # Calculated SNP weight (base_weight × zygosity_weight × snp_type_weight)
    "functional_modifier": 1.2,     # Modifier based on functional impact (e.g., from literature or annotation)
    "rarity_modifier": 1.1,         # 1.0 for common SNPs, 1.5 for rare, 2.0 for very rare
    "pathway_modifier": 1.2,        # Modifier for pathway importance (e.g., 1.2 if in a key pathway)
    "final_score": 0.4 * 1.2 * 1.0 * 1.1 * 1.2  # Product of all above modifiers
}

# SNP type weights
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

# Get SNP types using VEP
snp_results = get_snp_types(snp_list)

# --- Step 2:

# Check for invalid SNPs (not returned by VEP)
returned_snp_ids = {r['snp_id'] for r in snp_results}
input_snp_ids = set(snp_list.keys())
missing_snps = input_snp_ids - returned_snp_ids

if missing_snps:
    print("Invalid SNP supplied")
    exit(1)

# Extract unique gene symbols from SNPs (ignore 'UNKNOWN')
genes = list({r['gene'] for r in snp_results if r['gene'] != 'UNKNOWN'})

# Initialize STRING API client
string_api = StringAPI()

# Pull all downstream interaction genes for user's SNPs
invalid_genes = []
downstream_interactions = {}
for gene in genes:
    result = string_api.get_gene_interactions(gene, species=9606, required_score=400, include_visualization=False)
    interactions = result.get('interactions', [])
    if not interactions:
        invalid_genes.append(gene)
    else:
        downstream_interactions[gene] = interactions

if invalid_genes:
    print("Invalid gene supplied")
    exit(1)

#Print out the downstream interactions for each gene
for gene, interactions in downstream_interactions.items():
    print(f"\nDownstream interactions for {gene}:")
    for interaction in interactions:
        partner_id = interaction.get('stringId_B', 'UNKNOWN')
        partner_name = interaction.get('preferredName_B', partner_id)
        score = interaction.get('score', 0)
        if score >= 0.2:
            print(f"  Partner: {partner_name} (ID: {partner_id}), Confidence: {score:.2f}")

# --- Step 3:
# - Query for tissues affected by the SNPs



#
# ------------------------ Transcriptomic System -----------------------------------------------------------------------------
#   - mRNA expression levels
#   - Alternative splicing
#   - miRNA regulation
#   - note, differentially expressed genes
#   - Tools: GTEx, BioGPS, RNAseq datasets
#
# ------------------------ Proteomic System -----------------------------------------------------------------------------------
#   - Protein levels
#   - Protein-protein interactions
#   - Post-translational modifications (phosphorylation, glycosylation)
#   - Tools: STRING, BioGRID, ProteomicsDB
#
# ------------------------ Signaling and Receptor System -----------------------------------------------------------------------
#   - Receptors (e.g. α2A, D2, NMDA)
#   - Ligands (neurotransmitters, hormones)
#   - Intracellular signaling cascades
#   - Receptor desensitization/sensitization
#   - Tools: IUPHAR/BPS Guide to Pharmacology, Reactome, KEGG, NeuroMMSig
#
# ------------------------ Metabolic System -----------------------------------------------------------------------------------
#   - Enzymes (e.g. COMT, MAOA, MTHFR)
#   - Co-factors (e.g. B2, SAMe, copper)
#   - Pathways (methylation, detoxification, catecholamine degradation)
#   - Metabolites (dopamine, norepinephrine, GABA)
#   - Tools: KEGG, HMDB, WikiPathways
#
# ------------------------ Systems Biology Layers ------------------------------------------------------------------------------
#   - Organs and tissues (liver, brain, gut, etc.)
#   - Tissue-specific expression
#   - Cross-talk between systems (e.g. gut-brain axis, HPA axis)
#   - Transporters and barriers (e.g. blood-brain barrier, SLC transporters)
#   - Tools: Human Protein Atlas, Allen Brain Atlas, BRENDA
#
# Use network graphs linking SNP → gene → protein → pathway → phenotype
# Annotate receptor-level effects and neurotransmitter dynamics
# Track feedback loops and multisystem interdependencies

# Stage 3: Graph object with enriched nodes/edges 
#   !!!(this is currently just taking into account the SNPs and their downstream interactions, not the other systems)
#   - Node File (contains each node and all attributes)
#     - node_id	type	    layer	    location	expression	zygosity	receptor_type	neurotransmitter_affinity   tissues_relavence
#     - COMT	gene	    genomic	    PFC	        high	    true	    NA	            NA	                        ['example_tissue':1.0]
#     - α2A	    receptor	signaling	PFC	        medium	    false	    GPCR	        [norepinephrine	            ['example_tissue':0.5]
#   - Edge file
#     - source  target  type            confidence  directed
#     - COMT	DA      degradation     0.9         yes
#     - DA	    α2A	    activation	    0.6	        yes
#
# Pre-index node IDs
# Compress with .parquet or .feather if working in Python for speed
# Validate schema: no duplicate IDs, consistent types, edge directionality

# ------------------ MACHINE LEARNING PHASE ------------------


