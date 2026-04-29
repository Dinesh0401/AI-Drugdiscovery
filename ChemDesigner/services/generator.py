"""
Molecular generator service with simplified graph-based generation
"""
import random
import math
import logging
from typing import Dict, List, Any, Optional

class MolecularGenerator:
    """Simplified molecular generator for demonstration purposes"""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        
        # Predefined molecular fragments and scaffolds for realistic generation
        self.common_fragments = [
            'c1ccccc1',  # benzene
            'C1CCCCC1',  # cyclohexane
            'c1ccncc1',  # pyridine
            'c1ccc2ccccc2c1',  # naphthalene
            'C1CCC2CCCCC2C1',  # bicyclohexane
            'c1cnc2ccccc2c1',  # quinoline
            'c1ccc2[nH]ccc2c1',  # indole
            'C1CC2CCC1C2',  # norbornane
        ]
        
        self.functional_groups = [
            'O',     # hydroxyl
            'N',     # amino
            'C(=O)', # carbonyl
            'C(=O)O', # carboxyl
            'C#N',   # nitrile
            'S',     # thiol
            'C(F)(F)F', # trifluoromethyl
            'OC',    # methoxy
        ]
        
        # Drug-like molecules for sampling
        self.drug_like_molecules = [
            'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',  # Ibuprofen
            'CC(=O)OC1=CC=CC=C1C(=O)O',       # Aspirin
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # Caffeine
            'CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F', # Celecoxib
            'CC1=C(C(CCC1)(C)C)C=CC(=CC=CC(=CC(=O)O)C)C',  # Tretinoin
            'CN(C)CCCC1(C2=CC=CC=C2C3=CC=CC=C31)O',       # Dextromethorphan
        ]
    
    def generate_from_latent(self, latent_params: Dict[str, float]) -> Dict[str, Any]:
        """Generate molecule from latent space parameters"""
        try:
            # Use latent parameters to influence generation
            complexity = abs(latent_params.get('dim1', 0.0))
            polarity = latent_params.get('dim2', 0.0)
            ring_preference = latent_params.get('dim3', 0.0)
            heteroatom_preference = latent_params.get('dim4', 0.0)
            
            # Select base scaffold based on parameters
            scaffold = self._select_scaffold(complexity, ring_preference)
            
            # Add functional groups based on parameters
            molecule = self._modify_scaffold(scaffold, polarity, heteroatom_preference)
            
            return {
                'success': True,
                'smiles': molecule,
                'generation_params': latent_params,
                'structure_data': {
                    'scaffold_type': 'cyclic' if ring_preference > 0 else 'linear',
                    'complexity_level': min(10, int(complexity * 5)),
                    'polar_groups': int(abs(polarity) * 3)
                }
            }
            
        except Exception as e:
            self.logger.error(f"Error in latent generation: {e}")
            return {
                'success': False,
                'error': str(e)
            }
    
    def generate_random(self) -> Dict[str, Any]:
        """Generate a random drug-like molecule"""
        try:
            # Select random molecule from drug-like set
            base_molecule = random.choice(self.drug_like_molecules)
            
            # Optionally add small modifications
            if random.random() < 0.3:  # 30% chance of modification
                modified = self._add_simple_modification(base_molecule)
                if modified:
                    base_molecule = modified
            
            return {
                'success': True,
                'smiles': base_molecule
            }
            
        except Exception as e:
            self.logger.error(f"Error in random generation: {e}")
            return {
                'success': False,
                'error': str(e)
            }
    
    def optimize_for_properties(self, target_properties: Dict[str, float]) -> Dict[str, Any]:
        """Optimize molecule generation for target properties"""
        try:
            target_qed = target_properties.get('qed', 0.5)
            target_sa = target_properties.get('sa_score', 3.0)
            target_logp = target_properties.get('logp', 2.0)
            
            candidates = []
            
            # Generate multiple candidates and filter by properties
            for _ in range(50):
                # Bias generation parameters toward targets
                latent_params = self._generate_biased_latent(target_qed, target_sa, target_logp)
                result = self.generate_from_latent(latent_params)
                
                if result['success']:
                    candidates.append(result['smiles'])
            
            # Select best candidates (simplified selection)
            selected_molecules = list(set(candidates))[:10]  # Remove duplicates, limit to 10
            
            return {
                'success': True,
                'molecules': selected_molecules,
                'info': {
                    'total_generated': len(candidates),
                    'unique_molecules': len(selected_molecules),
                    'target_properties': target_properties
                }
            }
            
        except Exception as e:
            self.logger.error(f"Error in property optimization: {e}")
            return {
                'success': False,
                'error': str(e)
            }
    
    def _select_scaffold(self, complexity: float, ring_preference: float) -> str:
        """Select molecular scaffold based on parameters"""
        if ring_preference > 0.5:
            # Prefer cyclic structures
            if complexity > 0.7:
                return random.choice(self.common_fragments[3:])  # More complex rings
            else:
                return random.choice(self.common_fragments[:3])  # Simple rings
        else:
            # Linear or simpler structures
            if complexity > 0.5:
                return 'CCCCCC'  # Hexane chain
            else:
                return 'CCC'     # Propane
    
    def _modify_scaffold(self, scaffold: str, polarity: float, heteroatom_pref: float) -> str:
        """Add functional groups to scaffold"""
        molecule = scaffold
        
        # Add functional groups based on polarity
        if abs(polarity) > 0.3:
            num_groups = min(3, int(abs(polarity) * 4))
            for _ in range(num_groups):
                if heteroatom_pref > 0:
                    group = random.choice([g for g in self.functional_groups if 'N' in g or 'O' in g])
                else:
                    group = random.choice(self.functional_groups)
                
                # Simple attachment (not chemically rigorous)
                molecule = molecule + group
        
        return molecule
    
    def _add_simple_modification(self, smiles: str) -> Optional[str]:
        """Add simple modification to existing molecule"""
        modifications = ['C', 'O', 'F', 'Cl']
        
        try:
            modification = random.choice(modifications)
            # Simple concatenation (not chemically rigorous)
            return smiles + modification
        except:
            return None
    
    def _generate_biased_latent(self, target_qed: float, target_sa: float, target_logp: float) -> Dict[str, float]:
        """Generate latent parameters biased toward target properties"""
        
        # Rough mapping from properties to latent space (simplified)
        complexity_bias = max(0.1, min(0.9, 1.0 - target_qed))  # Lower QED = higher complexity
        polarity_bias = max(-0.5, min(0.5, (3.0 - target_logp) / 6.0))  # Lower LogP = more polar
        ring_bias = max(0.0, min(1.0, target_qed * 2 - 0.5))  # Higher QED = more rings
        heteroatom_bias = random.uniform(-0.3, 0.7)
        
        return {
            'dim1': complexity_bias + random.uniform(-0.2, 0.2),
            'dim2': polarity_bias + random.uniform(-0.1, 0.1),
            'dim3': ring_bias + random.uniform(-0.2, 0.2),
            'dim4': heteroatom_bias + random.uniform(-0.1, 0.1)
        }
