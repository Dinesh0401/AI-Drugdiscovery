"""
Property calculation service using RDKit for molecular descriptors
"""
import logging
import math
from typing import Dict, Any, Optional

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski, QED
    from rdkit.Contrib.SA_Score import sascorer
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available. Property calculations will use fallback methods.")

class PropertyCalculator:
    """Calculate molecular properties and drug-likeness metrics"""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        
        if not RDKIT_AVAILABLE:
            self.logger.warning("RDKit not available - using approximated calculations")
    
    def calculate_properties(self, smiles: str) -> Dict[str, Any]:
        """Calculate comprehensive molecular properties"""
        if not smiles or not isinstance(smiles, str):
            return self._get_empty_properties()
        
        try:
            if RDKIT_AVAILABLE:
                return self._calculate_with_rdkit(smiles)
            else:
                return self._calculate_fallback(smiles)
        except Exception as e:
            self.logger.error(f"Error calculating properties for {smiles}: {e}")
            return self._get_empty_properties()
    
    def _calculate_with_rdkit(self, smiles: str) -> Dict[str, Any]:
        """Calculate properties using RDKit"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return self._get_empty_properties()
        
        # Basic molecular properties
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        
        # Drug-likeness metrics
        qed_score = QED.qed(mol)
        
        # Synthetic accessibility score
        try:
            sa_score = sascorer.calculateScore(mol)
        except:
            sa_score = 5.0  # Default moderate difficulty
        
        # Lipinski's Rule of Five
        lipinski_violations = 0
        if mw > 500:
            lipinski_violations += 1
        if logp > 5:
            lipinski_violations += 1
        if hbd > 5:
            lipinski_violations += 1
        if hba > 10:
            lipinski_violations += 1
        
        # Additional descriptors
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        num_rings = Descriptors.RingCount(mol)
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        
        return {
            'molecular_weight': round(mw, 2),
            'logp': round(logp, 2),
            'hbd': hbd,
            'hba': hba,
            'tpsa': round(tpsa, 2),
            'rotatable_bonds': rotatable_bonds,
            'qed_score': round(qed_score, 3),
            'sa_score': round(sa_score, 2),
            'lipinski_violations': lipinski_violations,
            'num_atoms': num_atoms,
            'num_bonds': num_bonds,
            'num_rings': num_rings,
            'aromatic_rings': aromatic_rings,
            'drug_like': lipinski_violations <= 1 and qed_score >= 0.5,
            'synthetic_feasible': sa_score <= 6.0
        }
    
    def _calculate_fallback(self, smiles: str) -> Dict[str, Any]:
        """Fallback property calculation without RDKit"""
        # Very basic approximations based on SMILES string analysis
        
        # Count elements (very rough)
        carbon_count = smiles.count('C') + smiles.count('c')
        nitrogen_count = smiles.count('N') + smiles.count('n')
        oxygen_count = smiles.count('O') + smiles.count('o')
        sulfur_count = smiles.count('S') + smiles.count('s')
        
        # Rough molecular weight estimation
        mw = (carbon_count * 12.01 + nitrogen_count * 14.007 + 
              oxygen_count * 15.999 + sulfur_count * 32.065 + 
              smiles.count('H') * 1.008)
        
        # Very rough approximations
        logp = (carbon_count * 0.5 - nitrogen_count * 0.7 - oxygen_count * 1.2) / 2
        
        # Count potential H-bond donors/acceptors
        hbd = smiles.count('OH') + smiles.count('NH') + smiles.count('SH')
        hba = nitrogen_count + oxygen_count + sulfur_count
        
        # Rough estimates
        tpsa = (nitrogen_count * 12 + oxygen_count * 20)
        rotatable_bonds = smiles.count('-') // 2
        rings = smiles.count('1') + smiles.count('2') + smiles.count('3')
        
        # Simple drug-likeness approximation
        qed_score = max(0.0, min(1.0, 0.8 - abs(logp - 2) * 0.1 - max(0, mw - 500) * 0.001))
        sa_score = min(10.0, 3.0 + len(smiles) * 0.05)
        
        # Lipinski violations
        lipinski_violations = 0
        if mw > 500:
            lipinski_violations += 1
        if logp > 5:
            lipinski_violations += 1
        if hbd > 5:
            lipinski_violations += 1
        if hba > 10:
            lipinski_violations += 1
        
        return {
            'molecular_weight': round(mw, 2),
            'logp': round(logp, 2),
            'hbd': hbd,
            'hba': hba,
            'tpsa': round(tpsa, 2),
            'rotatable_bonds': rotatable_bonds,
            'qed_score': round(qed_score, 3),
            'sa_score': round(sa_score, 2),
            'lipinski_violations': lipinski_violations,
            'num_atoms': carbon_count + nitrogen_count + oxygen_count + sulfur_count,
            'num_bonds': max(0, len([c for c in smiles if c in 'CNOScnos']) - 1),
            'num_rings': rings,
            'aromatic_rings': min(rings, smiles.count('c')),
            'drug_like': lipinski_violations <= 1 and qed_score >= 0.5,
            'synthetic_feasible': sa_score <= 6.0
        }
    
    def _get_empty_properties(self) -> Dict[str, Any]:
        """Return empty property dictionary for invalid molecules"""
        return {
            'molecular_weight': 0.0,
            'logp': 0.0,
            'hbd': 0,
            'hba': 0,
            'tpsa': 0.0,
            'rotatable_bonds': 0,
            'qed_score': 0.0,
            'sa_score': 10.0,
            'lipinski_violations': 4,
            'num_atoms': 0,
            'num_bonds': 0,
            'num_rings': 0,
            'aromatic_rings': 0,
            'drug_like': False,
            'synthetic_feasible': False
        }
    
    def get_property_description(self, property_name: str) -> str:
        """Get human-readable description of a property"""
        descriptions = {
            'molecular_weight': 'Molecular weight in Daltons (Da)',
            'logp': 'Lipophilicity (LogP) - partition coefficient',
            'hbd': 'Hydrogen bond donors',
            'hba': 'Hydrogen bond acceptors',
            'tpsa': 'Topological polar surface area (Ų)',
            'rotatable_bonds': 'Number of rotatable bonds',
            'qed_score': 'Quantitative Estimate of Drug-likeness (0-1)',
            'sa_score': 'Synthetic Accessibility Score (1-10, lower is better)',
            'lipinski_violations': 'Lipinski Rule of Five violations',
            'num_atoms': 'Total number of atoms',
            'num_bonds': 'Total number of bonds',
            'num_rings': 'Number of rings',
            'aromatic_rings': 'Number of aromatic rings',
            'drug_like': 'Meets basic drug-likeness criteria',
            'synthetic_feasible': 'Synthetically feasible (SA ≤ 6)'
        }
        return descriptions.get(property_name, property_name)
