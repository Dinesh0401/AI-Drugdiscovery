/**
 * Molecular Discovery Platform - Main JavaScript Module
 * Handles molecular visualization, property display, and user interactions
 */

class MolecularPlatform {
    constructor() {
        this.molecules = new Map();
        this.currentMolecule = null;
        this.init();
    }
    
    init() {
        console.log('Molecular Discovery Platform initialized');
        this.bindEvents();
    }
    
    bindEvents() {
        // Navigation active state
        this.updateActiveNavigation();
        
        // Form validations
        this.setupFormValidations();
    }
    
    updateActiveNavigation() {
        const currentPath = window.location.pathname;
        const navLinks = document.querySelectorAll('.navbar-nav .nav-link');
        
        navLinks.forEach(link => {
            const href = link.getAttribute('href');
            if (href === currentPath) {
                link.classList.add('active');
            } else {
                link.classList.remove('active');
            }
        });
    }
    
    setupFormValidations() {
        const forms = document.querySelectorAll('form');
        forms.forEach(form => {
            form.addEventListener('submit', function(e) {
                if (!form.checkValidity()) {
                    e.preventDefault();
                    e.stopPropagation();
                }
                form.classList.add('was-validated');
            });
        });
    }
    
    // Utility functions for molecular operations
    static formatMolecularWeight(mw) {
        return `${parseFloat(mw).toFixed(2)} Da`;
    }
    
    static formatLogP(logp) {
        return parseFloat(logp).toFixed(2);
    }
    
    static formatPercentage(value) {
        return `${Math.round(value * 100)}%`;
    }
    
    static getDrugLikenessColor(qed) {
        if (qed >= 0.7) return 'success';
        if (qed >= 0.5) return 'warning';
        return 'danger';
    }
    
    static getSyntheticFeasibilityColor(sa) {
        if (sa <= 4) return 'success';
        if (sa <= 6) return 'warning';
        return 'danger';
    }
    
    // Animation utilities
    static animateValue(element, start, end, duration) {
        if (!element) return;
        
        const startTime = performance.now();
        const animate = (currentTime) => {
            const elapsed = currentTime - startTime;
            const progress = Math.min(elapsed / duration, 1);
            
            const current = start + (end - start) * progress;
            element.textContent = current.toFixed(2);
            
            if (progress < 1) {
                requestAnimationFrame(animate);
            }
        };
        
        requestAnimationFrame(animate);
    }
    
    static showToast(message, type = 'info') {
        const toast = document.createElement('div');
        toast.className = `toast align-items-center text-white bg-${type} border-0`;
        toast.setAttribute('role', 'alert');
        toast.innerHTML = `
            <div class="d-flex">
                <div class="toast-body">${message}</div>
                <button type="button" class="btn-close btn-close-white me-2 m-auto" data-bs-dismiss="toast"></button>
            </div>
        `;
        
        // Add to page
        let container = document.querySelector('.toast-container');
        if (!container) {
            container = document.createElement('div');
            container.className = 'toast-container position-fixed top-0 end-0 p-3';
            document.body.appendChild(container);
        }
        
        container.appendChild(toast);
        
        const bsToast = new bootstrap.Toast(toast);
        bsToast.show();
        
        // Auto-remove after shown
        toast.addEventListener('hidden.bs.toast', () => {
            toast.remove();
        });
    }
}

// Property calculation utilities
class PropertyCalculator {
    static calculateDrugLikeness(properties) {
        let score = 0;
        let maxScore = 5;
        
        // Molecular weight check
        if (properties.molecular_weight >= 150 && properties.molecular_weight <= 500) {
            score += 1;
        }
        
        // LogP check
        if (properties.logp >= -0.4 && properties.logp <= 5.6) {
            score += 1;
        }
        
        // Hydrogen bond checks
        if (properties.hbd <= 5) score += 1;
        if (properties.hba <= 10) score += 1;
        
        // TPSA check
        if (properties.tpsa <= 140) {
            score += 1;
        }
        
        return score / maxScore;
    }
    
    static getLipinskiViolations(properties) {
        let violations = 0;
        
        if (properties.molecular_weight > 500) violations++;
        if (properties.logp > 5) violations++;
        if (properties.hbd > 5) violations++;
        if (properties.hba > 10) violations++;
        
        return violations;
    }
    
    static getPropertyDescription(propertyName) {
        const descriptions = {
            'molecular_weight': 'The sum of atomic masses in the molecule',
            'logp': 'Measure of lipophilicity - how well the molecule dissolves in fats vs water',
            'hbd': 'Number of hydrogen bond donor groups (NH, OH)',
            'hba': 'Number of hydrogen bond acceptor atoms (N, O)',
            'tpsa': 'Polar surface area - important for membrane permeability',
            'rotatable_bonds': 'Number of bonds that can rotate freely - affects flexibility',
            'qed_score': 'Overall drug-likeness score from 0 (bad) to 1 (good)',
            'sa_score': 'Synthetic accessibility from 1 (easy) to 10 (very difficult)',
            'lipinski_violations': 'Number of Lipinski Rule of Five violations'
        };
        
        return descriptions[propertyName] || 'Property description not available';
    }
}

// SMILES utilities
class SMILESUtils {
    static isValid(smiles) {
        if (!smiles || typeof smiles !== 'string') return false;
        
        // Basic SMILES validation
        const invalidChars = /[<>|]/;
        if (invalidChars.test(smiles)) return false;
        
        // Check for balanced parentheses and brackets
        const brackets = smiles.match(/[\[\]]/g) || [];
        const parens = smiles.match(/[()]/g) || [];
        
        if (brackets.length % 2 !== 0 || parens.length % 2 !== 0) return false;
        
        return true;
    }
    
    static normalize(smiles) {
        return smiles.trim().replace(/\s+/g, '');
    }
    
    static getAtomCount(smiles) {
        // Simple approximation
        const atoms = smiles.match(/[A-Z][a-z]?/g) || [];
        return atoms.length;
    }
    
    static getBondCount(smiles) {
        // Simple approximation
        const bonds = smiles.match(/[-=#:]/g) || [];
        return bonds.length + Math.max(0, this.getAtomCount(smiles) - 1);
    }
}

// Initialize platform when DOM is loaded
document.addEventListener('DOMContentLoaded', function() {
    window.molecularPlatform = new MolecularPlatform();
});

// Export utilities for use in other modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        MolecularPlatform,
        PropertyCalculator,
        SMILESUtils
    };
}
