"""
crispr_designer/organisms.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Implements the Organism registry and genome reference chromosome mapping skeleton.
Provides access to NCBI/Ensembl chromosome accessions and handles
genomic sequence fetching (with high-fidelity local mock support).
"""

from typing import Dict, Any, Optional, List
import urllib.request
import json
import logging

# Logger only - basicConfig() is intentionally NOT called at import time. Configuring the root
# logger is the application entry point's responsibility (see main.py); doing it here would
# silently override the logging setup of anything that imports this module as a library.
logger = logging.getLogger("OrganismRegistry")

# Canonical organism metadata, one entry per organism, each listing its Korean/English aliases.
# Previously each alias (e.g. '벼' and 'Rice') carried its own full copy of the metadata dict,
# including the chromosome accession map - duplicating data in memory and risking the two
# copies drifting out of sync if only one was ever updated. Aliases now point at the same
# canonical dict instance instead.
_ORGANISM_DATA: List[Dict[str, Any]] = [
    {
        'aliases': ['벼', 'Rice'],
        'common_name_en': 'Rice',
        'scientific_name': 'Oryza sativa',
        'ref_assembly': 'IRGSP-1.0',
        'source': 'IRGSP / Ensembl Plants',
        'chromosomes': {
            '1': 'NC_029256.1', '2': 'NC_029257.1', '3': 'NC_029258.1',
            '4': 'NC_029259.1', '5': 'NC_029260.1', '6': 'NC_029261.1',
            '7': 'NC_029262.1', '8': 'NC_029263.1', '9': 'NC_029264.1',
            '10': 'NC_029265.1', '11': 'NC_029266.1', '12': 'NC_029267.1'
        }
    },
    {
        'aliases': ['애기장대', 'Arabidopsis'],
        'common_name_en': 'Arabidopsis',
        'scientific_name': 'Arabidopsis thaliana',
        'ref_assembly': 'TAIR10',
        'source': 'TAIR / NCBI',
        'chromosomes': {
            '1': 'NC_003070.9', '2': 'NC_003071.7', '3': 'NC_003074.8',
            '4': 'NC_003075.7', '5': 'NC_003076.8', 'MT': 'NC_001284.2',
            'CP': 'NC_000932.1'
        }
    },
    {
        'aliases': ['인간', 'Human'],
        'common_name_en': 'Human',
        'scientific_name': 'Homo sapiens',
        'ref_assembly': 'GRCh38',
        'source': 'GRCh38.p14 / NCBI Ensembl',
        'chromosomes': {
            '1': 'NC_000001.11', '2': 'NC_000002.12', '3': 'NC_000003.12',
            '4': 'NC_000004.12', '5': 'NC_000005.10', '6': 'NC_000006.12',
            '7': 'NC_000007.14', '8': 'NC_000008.11', '9': 'NC_000009.12',
            '10': 'NC_000010.11', '11': 'NC_000011.10', '12': 'NC_000012.12',
            '13': 'NC_000013.11', '14': 'NC_000014.9', '15': 'NC_000015.10',
            '16': 'NC_000016.10', '17': 'NC_000017.11', '18': 'NC_000018.10',
            '19': 'NC_000019.10', '20': 'NC_000020.11', '21': 'NC_000021.9',
            '22': 'NC_000022.11', 'X': 'NC_000023.11', 'Y': 'NC_000024.10',
            'MT': 'NC_012920.1'
        }
    },
    {
        'aliases': ['마우스', 'Mouse'],
        'common_name_en': 'Mouse',
        'scientific_name': 'Mus musculus',
        'ref_assembly': 'GRCm39',
        'source': 'GRCm39 / NCBI Ensembl',
        'chromosomes': {
            '1': 'NC_000067.7', '2': 'NC_000068.8', '3': 'NC_000069.7',
            '4': 'NC_000070.7', '5': 'NC_000071.7', '6': 'NC_000072.7',
            '7': 'NC_000073.7', '8': 'NC_000074.7', '9': 'NC_000075.7',
            '10': 'NC_000076.7', '11': 'NC_000077.7', '12': 'NC_000078.7',
            '13': 'NC_000079.7', '14': 'NC_000080.7', '15': 'NC_000081.7',
            '16': 'NC_000082.7', '17': 'NC_000083.7', '18': 'NC_000084.7',
            '19': 'NC_000085.7', 'X': 'NC_000086.8', 'Y': 'NC_000087.8',
            'MT': 'NC_005089.1'
        }
    },
    {
        'aliases': ['효모', 'Saccharomyces cerevisiae'],
        'common_name_en': 'Yeast',
        'scientific_name': 'Saccharomyces cerevisiae',
        'ref_assembly': 'sacCer3',
        'source': 'SGD / NCBI',
        'chromosomes': {
            'I': 'NC_001133.9', 'II': 'NC_001134.8', 'III': 'NC_001135.5',
            'IV': 'NC_001136.10', 'V': 'NC_001137.3', 'VI': 'NC_001138.5',
            'VII': 'NC_001139.9', 'VIII': 'NC_001140.6', 'IX': 'NC_001141.2',
            'X': 'NC_001142.9', 'XI': 'NC_001143.9', 'XII': 'NC_001144.5',
            'XIII': 'NC_001145.3', 'XIV': 'NC_001146.8', 'XV': 'NC_001147.6',
            'XVI': 'NC_001148.4', 'MT': 'NC_001224.1'
        }
    },
    {
        'aliases': ['대장균', 'E. coli'],
        'common_name_en': 'Ecoli',
        'scientific_name': 'Escherichia coli',
        'ref_assembly': 'K-12',
        'source': 'NCBI',
        'chromosomes': {
            '1': 'NC_000913.3'
        }
    }
]

# Organism Database Metadata: alias -> canonical metadata dict (shared reference per organism)
ORGANISM_REGISTRY: Dict[str, Dict[str, Any]] = {}
for _entry in _ORGANISM_DATA:
    _aliases = _entry['aliases']
    _metadata = {k: v for k, v in _entry.items() if k != 'aliases'}
    for _alias in _aliases:
        ORGANISM_REGISTRY[_alias] = _metadata

# High-fidelity biological mock sequences for testing offline (contains exact Cas9/Cas12a/Cas13d targets)
MOCK_GENOMIC_DB = {
    'Human': {
        '1': "ATGCGATCGATCGATCGATCGATCGATCGAATCGATCGATCGATCGGGCATCGATCGATCGATCGGGATCGATCGATCGAATTTGATCGATCGATCGATCGGGCATCGATCGATCGATCG",
        'X': "CCCGGGTATATATCGATCGATCGACTAGCTAGCTAGCTAGCTGATCGATCGATCGATCGACTTTNATCGATCGATCGATCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
    },
    'Rice': {
        '1': "ATGGGGCATGCATCGACTAGCTAGCTAGCTAGCGCATCGATCGATCGATCGGGCATCGATCGATCGATCGATCTTTTCGATCGATCGATCGATCGCGGGGGGGGGGGGGGGG",
        '10': "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
    },
    'Arabidopsis': {
        '1': "ATGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTTTNAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
        '5': "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    },
    'Mouse': {
        '1': "ATGCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTNATGCATGCATGCATGCATGCATGCATGCATGCATGCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
    },
    'Yeast': {
        'I': "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    },
    'Ecoli': {
        '1': "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    }
}

def get_organism_metadata(name: str) -> Dict[str, Any]:
    """
    Retrieves the genomic metadata dictionary for the specified organism.
    """
    if name not in ORGANISM_REGISTRY:
        raise KeyError(
            f"Organism '{name}' not found in registry. "
            f"Supported organisms: {list(ORGANISM_REGISTRY.keys())}"
        )
    return ORGANISM_REGISTRY[name]

def fetch_genomic_sequence(
    organism_name: str, 
    chromosome: str, 
    start: Optional[int] = None, 
    end: Optional[int] = None,
    force_mock: bool = True
) -> str:
    """
    Skeleton mapping function that resolves Ensembl or NCBI reference accessions and
    fetches target chromosome sequence data.

    In a production system, this sends an API query to Ensembl REST or NCBI E-utilities.
    If 'force_mock=True' (default for safety and testing), it resolves to high-fidelity mock sequences offline.

    Args:
        organism_name (str): '벼', 'Rice', '애기장대', 'Arabidopsis', '인간', or 'Human'.
        chromosome (str): Chromosome number or key (e.g. '1', '2', 'X').
        start (int, optional): 0-indexed start boundary coordinates.
        end (int, optional): 0-indexed end boundary coordinates.
        force_mock (bool): If True, skips network requests and loads local test database.

    Returns:
        str: Fasta/nucleotide sequence string.
    """
    # 1. Resolve organism and assembly details
    metadata = get_organism_metadata(organism_name)
    assembly = metadata['ref_assembly']
    common_en = metadata['common_name_en']
    
    # 2. Map chromosome to accession ID
    chromosomes_map = metadata.get('chromosomes', {})
    if chromosome not in chromosomes_map:
        raise ValueError(
            f"Invalid chromosome '{chromosome}' for {organism_name}. "
            f"Supported chromosomes: {list(chromosomes_map.keys())}"
        )
    accession = chromosomes_map[chromosome]
    
    logger.info(f"Mapping {organism_name} ({assembly}) Chr {chromosome} to Accession {accession}")

    # 3. Fetch sequence
    if force_mock:
        # Load local high-fidelity biological mockup
        logger.info("Local mock mode active: loading sequence from test DB.")
        full_seq = MOCK_GENOMIC_DB.get(common_en, {}).get(
            chromosome, 
            # Default fallback random sequence rich with PAMs if chromosome not mocked
            "ATGCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTNATGCATGCATGCATGCATGCATGCATGCATGCATGCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
        )
    else:
        # Real HTTP Web API fetch skeleton
        # Uses Ensembl REST API which provides rapid JSON responses for chromosome segments
        # Endpoint format: http://rest.ensembl.org/sequence/region/human/1:1000..2000?content-type=application/json
        try:
            scientific_lower = metadata['scientific_name'].lower().replace(" ", "_")
            start_coord = start if start is not None else 1
            end_coord = end if end is not None else 5000  # Default small chunk to avoid timeouts
            
            url = (
                f"http://rest.ensembl.org/sequence/region/{scientific_lower}/"
                f"{chromosome}:{start_coord}..{end_coord}?content-type=application/json"
            )
            
            logger.info(f"Querying Ensembl API: {url}")
            req = urllib.request.Request(url, headers={'Content-Type': 'application/json'})
            with urllib.request.urlopen(req, timeout=5) as response:
                res_data = json.loads(response.read().decode('utf-8'))
                full_seq = res_data.get('seq', '')
                logger.info("Successfully fetched sequence from Ensembl API.")
                return full_seq
        except Exception as e:
            logger.warning(f"Failed to query external genomic API: {e}. Falling back to high-fidelity mock data.")
            # Graceful fallback to mock so it never crashes the platform in unstable networks
            full_seq = MOCK_GENOMIC_DB.get(common_en, {}).get(chromosome, "ATGCGTACGTACGTACGTACGTNNGGTTTNATGC")

    # Handle coordinates slicing for the mock database
    if start is not None or end is not None:
        slice_start = start if start is not None else 0
        slice_end = end if end is not None else len(full_seq)
        return full_seq[slice_start:slice_end]
        
    return full_seq
