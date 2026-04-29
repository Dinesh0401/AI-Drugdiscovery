/**
 * Canvas-based molecular structure renderer
 * Provides 2D molecular visualization on HTML5 Canvas
 */

class MoleculeRenderer {
    constructor(canvas) {
        this.canvas = canvas;
        this.ctx = canvas.getContext('2d');
        this.width = canvas.width;
        this.height = canvas.height;
        
        // Rendering parameters
        this.atomRadius = 12;
        this.bondWidth = 3;
        this.colors = {
            carbon: '#666666',
            nitrogen: '#3050f8',
            oxygen: '#ff0d0d',
            sulfur: '#ffff30',
            phosphorus: '#ff8000',
            halogen: '#1dc51d',
            hydrogen: '#ffffff',
            other: '#ff1dce',
            bond: '#cccccc',
            background: '#1a1a1a'
        };
        
        this.currentMolecule = null;
        this.setupCanvas();
    }
    
    setupCanvas() {
        // Set up canvas with high DPI support
        const devicePixelRatio = window.devicePixelRatio || 1;
        const rect = this.canvas.getBoundingClientRect();
        
        this.canvas.width = rect.width * devicePixelRatio;
        this.canvas.height = rect.height * devicePixelRatio;
        
        this.ctx.scale(devicePixelRatio, devicePixelRatio);
        this.canvas.style.width = rect.width + 'px';
        this.canvas.style.height = rect.height + 'px';
        
        this.width = rect.width;
        this.height = rect.height;
        
        this.clearCanvas();
    }
    
    clearCanvas() {
        this.ctx.fillStyle = this.colors.background;
        this.ctx.fillRect(0, 0, this.width, this.height);
    }
    
    renderMolecule(smiles, properties = null) {
        this.currentMolecule = { smiles, properties };
        this.clearCanvas();
        
        if (!smiles) {
            this.renderPlaceholder();
            return;
        }
        
        try {
            // Generate simplified 2D coordinates
            const structure = this.generateStructure(smiles, properties);
            this.drawMolecule(structure);
        } catch (error) {
            console.warn('Error rendering molecule:', error);
            this.renderError(smiles);
        }
    }
    
    generateStructure(smiles, properties) {
        // Simplified structure generation based on SMILES and properties
        const atoms = this.parseAtoms(smiles);
        const bonds = this.parseBonds(smiles);
        
        // Generate 2D layout
        const coordinates = this.generate2DLayout(atoms, bonds, properties);
        
        return {
            atoms: atoms,
            bonds: bonds,
            coordinates: coordinates
        };
    }
    
    parseAtoms(smiles) {
        const atoms = [];
        const atomRegex = /[A-Z][a-z]?/g;
        let match;
        
        while ((match = atomRegex.exec(smiles)) !== null) {
            atoms.push({
                symbol: match[0],
                index: atoms.length,
                position: match.index
            });
        }
        
        // Add implicit hydrogens for visualization
        if (atoms.length === 0) {
            atoms.push({ symbol: 'C', index: 0, position: 0 });
        }
        
        return atoms;
    }
    
    parseBonds(smiles) {
        const bonds = [];
        
        // Simple bond detection (very basic)
        const bondChars = smiles.match(/[-=#]/g) || [];
        
        // Create linear bonds for simple visualization
        for (let i = 0; i < Math.max(1, bondChars.length); i++) {
            bonds.push({
                from: i,
                to: i + 1,
                type: bondChars[i] || 'single',
                aromatic: smiles.includes('c') && i < 6 // Simple aromatic detection
            });
        }
        
        return bonds;
    }
    
    generate2DLayout(atoms, bonds, properties) {
        const coordinates = [];
        const centerX = this.width / 2;
        const centerY = this.height / 2;
        
        if (atoms.length === 1) {
            // Single atom
            coordinates.push({ x: centerX, y: centerY });
            return coordinates;
        }
        
        // Determine layout based on properties
        const isRingSystem = properties && (properties.num_rings > 0 || properties.aromatic_rings > 0);
        const isLinear = !isRingSystem && atoms.length <= 6;
        
        if (isRingSystem) {
            // Ring layout
            this.generateRingLayout(atoms, coordinates, centerX, centerY);
        } else if (isLinear) {
            // Linear layout
            this.generateLinearLayout(atoms, coordinates, centerX, centerY);
        } else {
            // Complex branched layout
            this.generateBranchedLayout(atoms, coordinates, centerX, centerY);
        }
        
        return coordinates;
    }
    
    generateRingLayout(atoms, coordinates, centerX, centerY) {
        const radius = Math.min(this.width, this.height) * 0.2;
        const angleStep = (2 * Math.PI) / atoms.length;
        
        for (let i = 0; i < atoms.length; i++) {
            const angle = i * angleStep - Math.PI / 2; // Start from top
            const x = centerX + radius * Math.cos(angle);
            const y = centerY + radius * Math.sin(angle);
            coordinates.push({ x, y });
        }
    }
    
    generateLinearLayout(atoms, coordinates, centerX, centerY) {
        const totalWidth = Math.min(this.width * 0.6, atoms.length * 40);
        const startX = centerX - totalWidth / 2;
        const stepX = totalWidth / Math.max(1, atoms.length - 1);
        
        for (let i = 0; i < atoms.length; i++) {
            const x = startX + i * stepX;
            const y = centerY + (i % 2 === 1 ? 20 : 0); // Slight zigzag for visual interest
            coordinates.push({ x, y });
        }
    }
    
    generateBranchedLayout(atoms, coordinates, centerX, centerY) {
        // Grid-based layout for complex molecules
        const cols = Math.ceil(Math.sqrt(atoms.length));
        const rows = Math.ceil(atoms.length / cols);
        const spacing = 50;
        
        const totalWidth = (cols - 1) * spacing;
        const totalHeight = (rows - 1) * spacing;
        const startX = centerX - totalWidth / 2;
        const startY = centerY - totalHeight / 2;
        
        for (let i = 0; i < atoms.length; i++) {
            const row = Math.floor(i / cols);
            const col = i % cols;
            
            const x = startX + col * spacing;
            const y = startY + row * spacing;
            coordinates.push({ x, y });
        }
    }
    
    drawMolecule(structure) {
        // Draw bonds first (so they appear behind atoms)
        this.drawBonds(structure);
        
        // Draw atoms
        this.drawAtoms(structure);
        
        // Draw molecule info
        this.drawMoleculeInfo(structure);
    }
    
    drawBonds(structure) {
        const { bonds, coordinates } = structure;
        
        this.ctx.strokeStyle = this.colors.bond;
        this.ctx.lineWidth = this.bondWidth;
        this.ctx.lineCap = 'round';
        
        for (const bond of bonds) {
            if (bond.from >= coordinates.length || bond.to >= coordinates.length) continue;
            
            const fromCoord = coordinates[bond.from];
            const toCoord = coordinates[bond.to];
            
            if (!fromCoord || !toCoord) continue;
            
            this.ctx.beginPath();
            this.ctx.moveTo(fromCoord.x, fromCoord.y);
            this.ctx.lineTo(toCoord.x, toCoord.y);
            this.ctx.stroke();
            
            // Draw double/triple bonds
            if (bond.type === '=' || bond.type === '#') {
                this.drawMultipleBond(fromCoord, toCoord, bond.type);
            }
        }
    }
    
    drawMultipleBond(from, to, bondType) {
        const offset = bondType === '=' ? 4 : 6;
        const dx = to.x - from.x;
        const dy = to.y - from.y;
        const length = Math.sqrt(dx * dx + dy * dy);
        
        if (length === 0) return;
        
        const perpX = -dy / length * offset;
        const perpY = dx / length * offset;
        
        // Draw additional bond line(s)
        this.ctx.beginPath();
        this.ctx.moveTo(from.x + perpX, from.y + perpY);
        this.ctx.lineTo(to.x + perpX, to.y + perpY);
        this.ctx.stroke();
        
        if (bondType === '#') {
            this.ctx.beginPath();
            this.ctx.moveTo(from.x - perpX, from.y - perpY);
            this.ctx.lineTo(to.x - perpX, to.y - perpY);
            this.ctx.stroke();
        }
    }
    
    drawAtoms(structure) {
        const { atoms, coordinates } = structure;
        
        for (let i = 0; i < atoms.length && i < coordinates.length; i++) {
            const atom = atoms[i];
            const coord = coordinates[i];
            
            if (!coord) continue;
            
            // Draw atom circle
            const color = this.getAtomColor(atom.symbol);
            this.ctx.fillStyle = color;
            this.ctx.beginPath();
            this.ctx.arc(coord.x, coord.y, this.atomRadius, 0, 2 * Math.PI);
            this.ctx.fill();
            
            // Draw atom symbol
            this.ctx.fillStyle = atom.symbol === 'C' ? '#ffffff' : '#000000';
            this.ctx.font = 'bold 12px Arial';
            this.ctx.textAlign = 'center';
            this.ctx.textBaseline = 'middle';
            this.ctx.fillText(atom.symbol, coord.x, coord.y);
        }
    }
    
    drawMoleculeInfo(structure) {
        if (!this.currentMolecule || !this.currentMolecule.properties) return;
        
        const properties = this.currentMolecule.properties;
        
        // Draw property info in corner
        this.ctx.fillStyle = 'rgba(0, 0, 0, 0.7)';
        this.ctx.fillRect(10, 10, 160, 80);
        
        this.ctx.fillStyle = '#ffffff';
        this.ctx.font = '11px Arial';
        this.ctx.textAlign = 'left';
        this.ctx.textBaseline = 'top';
        
        const info = [
            `Atoms: ${properties.num_atoms || structure.atoms.length}`,
            `Bonds: ${properties.num_bonds || structure.bonds.length}`,
            `MW: ${properties.molecular_weight || 'N/A'} Da`,
            `LogP: ${properties.logp || 'N/A'}`,
            `QED: ${properties.qed_score || 'N/A'}`
        ];
        
        info.forEach((text, i) => {
            this.ctx.fillText(text, 15, 15 + i * 14);
        });
    }
    
    getAtomColor(symbol) {
        switch (symbol.toLowerCase()) {
            case 'c': return this.colors.carbon;
            case 'n': return this.colors.nitrogen;
            case 'o': return this.colors.oxygen;
            case 's': return this.colors.sulfur;
            case 'p': return this.colors.phosphorus;
            case 'f': case 'cl': case 'br': case 'i':
                return this.colors.halogen;
            case 'h': return this.colors.hydrogen;
            default: return this.colors.other;
        }
    }
    
    renderPlaceholder() {
        this.ctx.fillStyle = '#666666';
        this.ctx.font = '16px Arial';
        this.ctx.textAlign = 'center';
        this.ctx.textBaseline = 'middle';
        this.ctx.fillText('Generate a molecule to see structure', this.width / 2, this.height / 2);
        
        // Draw placeholder molecular icon
        this.drawPlaceholderMolecule();
    }
    
    drawPlaceholderMolecule() {
        const centerX = this.width / 2;
        const centerY = this.height / 2 - 40;
        const radius = 8;
        
        // Draw simple benzene ring as placeholder
        const atoms = 6;
        const ringRadius = 30;
        
        // Draw bonds
        this.ctx.strokeStyle = '#444444';
        this.ctx.lineWidth = 2;
        
        for (let i = 0; i < atoms; i++) {
            const angle1 = (i * 2 * Math.PI) / atoms - Math.PI / 2;
            const angle2 = ((i + 1) % atoms * 2 * Math.PI) / atoms - Math.PI / 2;
            
            const x1 = centerX + ringRadius * Math.cos(angle1);
            const y1 = centerY + ringRadius * Math.sin(angle1);
            const x2 = centerX + ringRadius * Math.cos(angle2);
            const y2 = centerY + ringRadius * Math.sin(angle2);
            
            this.ctx.beginPath();
            this.ctx.moveTo(x1, y1);
            this.ctx.lineTo(x2, y2);
            this.ctx.stroke();
        }
        
        // Draw atoms
        this.ctx.fillStyle = this.colors.carbon;
        for (let i = 0; i < atoms; i++) {
            const angle = (i * 2 * Math.PI) / atoms - Math.PI / 2;
            const x = centerX + ringRadius * Math.cos(angle);
            const y = centerY + ringRadius * Math.sin(angle);
            
            this.ctx.beginPath();
            this.ctx.arc(x, y, radius, 0, 2 * Math.PI);
            this.ctx.fill();
        }
    }
    
    renderError(smiles) {
        this.ctx.fillStyle = '#ff4444';
        this.ctx.font = '14px Arial';
        this.ctx.textAlign = 'center';
        this.ctx.textBaseline = 'middle';
        this.ctx.fillText('Error rendering molecule', this.width / 2, this.height / 2);
        
        this.ctx.fillStyle = '#888888';
        this.ctx.font = '12px Arial';
        this.ctx.fillText(`SMILES: ${smiles.substring(0, 30)}...`, this.width / 2, this.height / 2 + 25);
    }
    
    // Public methods for interaction
    exportCanvas() {
        return this.canvas.toDataURL('image/png');
    }
    
    resize() {
        this.setupCanvas();
        if (this.currentMolecule) {
            this.renderMolecule(this.currentMolecule.smiles, this.currentMolecule.properties);
        }
    }
}

// Handle canvas resize on window resize
window.addEventListener('resize', function() {
    const canvases = document.querySelectorAll('canvas.molecule-canvas');
    canvases.forEach(canvas => {
        if (canvas.moleculeRenderer) {
            canvas.moleculeRenderer.resize();
        }
    });
});

// Export for use in other modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = MoleculeRenderer;
}
