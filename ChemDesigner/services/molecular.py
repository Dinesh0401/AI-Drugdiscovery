"""
Molecular service for handling SMILES validation and basic molecular operations
"""
import logging
from typing import Optional, Dict, Any

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available. Some functionality will be limited.")

class MolecularService:
    """Service for basic molecular operations"""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        
        # Common chemical formulas to SMILES mapping for user convenience
        self.formula_to_smiles = {
            'H2O': 'O',
            'CO2': 'O=C=O', 
            'CO': '[C-]#[O+]',
            'O2': 'O=O',
            'N2': 'N#N',
            'H2': '[H][H]',
            'NH3': 'N',
            'CH4': 'C',
            'C2H6': 'CC',
            'C2H4': 'C=C',
            'C2H2': 'C#C',
            'C6H6': 'c1ccccc1',  # Benzene
            'C6H12O6': 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',  # Glucose
            'C8H10N4O2': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # Caffeine
            'C9H8O4': 'CC(=O)OC1=CC=CC=C1C(=O)O',  # Aspirin
            'NaCl': '[Na+].[Cl-]',
            'CaCO3': '[Ca+2].[O-]C([O-])=O',
            'H2SO4': 'OS(=O)(=O)O'
        }
        
        if not RDKIT_AVAILABLE:
            self.logger.warning("RDKit not available - using fallback implementations")
    
    def normalize_input(self, user_input: str) -> tuple[str, str]:
        """
        Normalize user input to SMILES format
        Returns: (smiles, input_type) where input_type is 'formula', 'smiles', or 'unknown'
        """
        if not user_input or not isinstance(user_input, str):
            return "", "unknown"
        
        # Clean the input
        clean_input = user_input.strip().replace(' ', '')
        
        # Handle Unicode subscripts (CO₂ -> CO2)
        unicode_map = {'₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4', 
                      '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9'}
        for unicode_char, normal_char in unicode_map.items():
            clean_input = clean_input.replace(unicode_char, normal_char)
        
        # Check if it's a known formula
        if clean_input.upper() in self.formula_to_smiles:
            return self.formula_to_smiles[clean_input.upper()], "formula"
        
        # Try common case variations
        for formula, smiles in self.formula_to_smiles.items():
            if clean_input.upper() == formula.upper():
                return smiles, "formula"
        
        # If it looks like a SMILES string (contains lowercase letters, brackets, etc.)
        if any(c in clean_input for c in ['c', '[', ']', '(', ')', '=', '#', '@']):
            return clean_input, "smiles"
        
        # Check if it might be a simple molecular formula pattern
        import re
        if re.match(r'^[A-Z][a-z]?(\d*[A-Z][a-z]?\d*)*$', clean_input):
            # It looks like a molecular formula but we don't recognize it
            return clean_input, "unknown_formula"
        
        # Assume it's a SMILES string
        return clean_input, "smiles"
    
    def is_valid_smiles(self, smiles: str) -> bool:
        """Validate a SMILES string"""
        if not smiles or not isinstance(smiles, str):
            return False
            
        if not RDKIT_AVAILABLE:
            # Basic validation without RDKit
            return len(smiles.strip()) > 0 and not any(char in smiles for char in ['<', '>', '|'])
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception as e:
            self.logger.debug(f"SMILES validation failed: {e}")
            return False
    
    def canonicalize_smiles(self, smiles: str) -> Optional[str]:
        """Convert SMILES to canonical form"""
        if not self.is_valid_smiles(smiles):
            return None
            
        if not RDKIT_AVAILABLE:
            return smiles.strip()
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return Chem.MolToSmiles(mol, canonical=True)
        except Exception as e:
            self.logger.error(f"Error canonicalizing SMILES: {e}")
            return None
    
    def get_molecular_formula(self, smiles: str) -> Optional[str]:
        """Get molecular formula from SMILES"""
        if not self.is_valid_smiles(smiles):
            return None
            
        if not RDKIT_AVAILABLE:
            return "C?H?N?O?"  # Placeholder when RDKit unavailable
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return Chem.rdMolDescriptors.CalcMolFormula(mol)
        except Exception as e:
            self.logger.error(f"Error getting molecular formula: {e}")
            return None
    
    def get_atom_count(self, smiles: str) -> int:
        """Get total atom count"""
        if not self.is_valid_smiles(smiles):
            return 0
            
        if not RDKIT_AVAILABLE:
            # Simple approximation
            return len([c for c in smiles if c.isupper()])
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return 0
            return mol.GetNumAtoms()
        except Exception as e:
            self.logger.error(f"Error getting atom count: {e}")
            return 0
    
    def get_bond_count(self, smiles: str) -> int:
        """Get total bond count"""
        if not self.is_valid_smiles(smiles):
            return 0
            
        if not RDKIT_AVAILABLE:
            # Simple approximation
            return max(0, len([c for c in smiles if c.isupper()]) - 1)
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return 0
            return mol.GetNumBonds()
        except Exception as e:
            self.logger.error(f"Error getting bond count: {e}")
            return 0
    
    def smiles_to_graph_data(self, smiles: str) -> Dict[str, Any]:
        """Convert SMILES to simplified graph representation"""
        if not self.is_valid_smiles(smiles):
            return {}
        
        try:
            if RDKIT_AVAILABLE:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    return {}
                
                # Extract basic graph information
                atoms = []
                bonds = []
                
                for atom in mol.GetAtoms():
                    atoms.append({
                        'symbol': atom.GetSymbol(),
                        'atomic_num': atom.GetAtomicNum(),
                        'formal_charge': atom.GetFormalCharge(),
                        'hybridization': str(atom.GetHybridization()),
                        'index': atom.GetIdx()
                    })
                
                for bond in mol.GetBonds():
                    bonds.append({
                        'begin_atom': bond.GetBeginAtomIdx(),
                        'end_atom': bond.GetEndAtomIdx(),
                        'bond_type': str(bond.GetBondType()),
                        'is_aromatic': bond.GetIsAromatic()
                    })
                
                return {
                    'atoms': atoms,
                    'bonds': bonds,
                    'num_atoms': len(atoms),
                    'num_bonds': len(bonds)
                }
            else:
                # Fallback without RDKit
                return {
                    'atoms': [],
                    'bonds': [],
                    'num_atoms': self.get_atom_count(smiles),
                    'num_bonds': self.get_bond_count(smiles)
                }
                
        except Exception as e:
            self.logger.error(f"Error converting SMILES to graph: {e}")
            return {}
