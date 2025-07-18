import requests
import json
from typing import List, Dict, Optional, Union, Tuple

class StringAPI:
    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize the StringAPI client.
        
        Args:
            api_key (str, optional): Your STRING API key. If not provided, will use public access.
        """
        self.base_url = "https://string-db.org/api"
        self.api_key = api_key

    def get_gene_interactions(self, gene: str, species: int = 9606,
                            required_score: int = 400,
                            include_visualization: bool = True) -> Dict[str, Union[Dict, str]]:
        """
        Get all direct interactions for a single gene, including optional visualization.
        
        Args:
            gene (str): Gene identifier (e.g., gene name, Ensembl ID)
            species (int): NCBI taxonomy identifier (default: 9606 for human)
            required_score (int): Minimum interaction score (0-1000, default: 400)
            include_visualization (bool): Whether to include network visualization URL (default: True)
            
        Returns:
            Dict containing:
                - 'interactions': Dict of interaction data
                - 'visualization_url': str URL to network visualization (if include_visualization is True)
        """
        # Get interaction partners
        endpoint = f"{self.base_url}/json/interaction_partners"
        
        params = {
            "identifier": gene,
            "species": species,
            "required_score": required_score
        }
        
        if self.api_key:
            params["api_key"] = self.api_key
            
        response = requests.get(endpoint, params=params)
        response.raise_for_status()
        interactions = response.json()
        
        result = {
            'interactions': interactions
        }
        
        # Get visualization URL if requested
        if include_visualization:
            # Create list of genes including the query gene and its partners
            all_genes = [gene] + [interaction['stringId_B'] for interaction in interactions]
            result['visualization_url'] = self.get_network_image_url(all_genes, species, required_score)
        
        return result

    def get_network_image_url(self, genes: List[str], species: int = 9606,
                            required_score: int = 400, network_type: str = "functional",
                            additional_network_nodes: int = 0,
                            image_format: str = "image") -> str:
        """
        Get the URL for the visual network image of gene interactions.
        
        Args:
            genes (List[str]): List of gene identifiers
            species (int): NCBI taxonomy identifier (default: 9606 for human)
            required_score (int): Minimum interaction score (0-1000, default: 400)
            network_type (str): Type of network ("functional" or "physical", default: "functional")
            additional_network_nodes (int): Number of additional nodes to include (default: 0)
            image_format (str): Format of the image ("image" or "highres_image", default: "image")
            
        Returns:
            str: URL to the network image
        """
        endpoint = f"{self.base_url}/{image_format}/network"
        
        params = {
            "identifiers": "%0d".join(genes),
            "species": species,
            "required_score": required_score,
            "network_type": network_type,
            "additional_network_nodes": additional_network_nodes
        }
        
        if self.api_key:
            params["api_key"] = self.api_key
            
        # Build the URL with parameters
        query_string = "&".join(f"{k}={v}" for k, v in params.items())
        return f"{endpoint}?{query_string}"

    def get_downstream_interactions_for_genes(self, genes: List[str], species: int = 9606, required_score: int = 400) -> Dict[str, List[Dict]]:
        """
        For a list of genes, get all direct downstream interaction partners and their confidence scores.
        Returns a dict mapping gene -> list of interaction dicts.
        """
        downstream = {}
        for gene in genes:
            try:
                result = self.get_gene_interactions(gene, species=species, required_score=required_score, include_visualization=False)
                downstream[gene] = result['interactions']
            except Exception as e:
                downstream[gene] = []
        return downstream

    @staticmethod
    def extract_interaction_partners(interactions: List[Dict], min_confidence: float = 0.2) -> List[Tuple[str, float]]:
        """
        Given a list of STRING interaction dicts, extract (partner_gene, confidence) tuples above min_confidence (0-1 scale).
        """
        partners = []
        for interaction in interactions:
            partner = interaction.get('stringId_B')
            score = interaction.get('score', 0)
            if score >= min_confidence:
                partners.append((partner, score))
        return partners

# Example usage
if __name__ == "__main__":
    # Initialize the API client (with or without API key)
    string_api = StringAPI()
    
    # Example: Get all interactions for a single gene
    gene = "MAOA"
    result = string_api.get_gene_interactions(gene)
    
    # Print interaction data
    print(f"\nDirect interactions for {gene}:")
    print(json.dumps(result['interactions'], indent=2))
    
    # Print visualization URL
    print(f"\nNetwork visualization URL:")
    print(result['visualization_url']) 