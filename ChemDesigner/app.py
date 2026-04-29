import os
import logging
from flask import Flask, render_template, request, jsonify, flash, redirect, url_for
from services.molecular import MolecularService
from services.properties import PropertyCalculator
from services.generator import MolecularGenerator

# Configure logging
logging.basicConfig(level=logging.DEBUG)

# Create Flask app
app = Flask(__name__)
app.secret_key = os.environ.get("SESSION_SECRET", "molecular-discovery-key-2024")

# Initialize services
molecular_service = MolecularService()
property_calculator = PropertyCalculator()
molecular_generator = MolecularGenerator()

@app.route('/')
def index():
    """Home page with overview of molecular generation system"""
    return render_template('index.html')

@app.route('/generate')
def generate_page():
    """Interactive molecular generation page"""
    return render_template('generate.html')

@app.route('/analyze')
def analyze_page():
    """Molecular analysis and property calculation page"""
    return render_template('analyze.html')

@app.route('/api/generate_molecule', methods=['POST'])
def generate_molecule():
    """Generate a molecule based on latent space parameters"""
    try:
        data = request.get_json() if request.is_json else request.form
        
        # Extract parameters
        latent_params = {
            'dim1': float(data.get('dim1', 0.0)),
            'dim2': float(data.get('dim2', 0.0)),
            'dim3': float(data.get('dim3', 0.0)),
            'dim4': float(data.get('dim4', 0.0))
        }
        
        # Generate molecule
        result = molecular_generator.generate_from_latent(latent_params)
        
        if result['success']:
            # Calculate properties
            smiles = result['smiles']
            properties = property_calculator.calculate_properties(smiles)
            
            return jsonify({
                'success': True,
                'smiles': smiles,
                'properties': properties,
                'structure_data': result.get('structure_data', {})
            })
        else:
            return jsonify({
                'success': False,
                'error': result.get('error', 'Generation failed')
            })
            
    except Exception as e:
        app.logger.error(f"Error in generate_molecule: {str(e)}")
        return jsonify({'success': False, 'error': str(e)})

@app.route('/api/calculate_properties', methods=['POST'])
def calculate_properties():
    """Calculate properties for a given SMILES string or chemical formula"""
    try:
        data = request.get_json() if request.is_json else request.form
        user_input = data.get('smiles', '').strip()
        
        if not user_input:
            return jsonify({'success': False, 'error': 'Chemical formula or SMILES string required'})
        
        # Normalize input (convert formulas like CO2 to SMILES)
        smiles, input_type = molecular_service.normalize_input(user_input)
        
        if input_type == "unknown_formula":
            return jsonify({
                'success': False, 
                'error': f'Unknown chemical formula: {user_input}. Try entering a SMILES string instead.'
            })
        
        # Validate SMILES
        if not molecular_service.is_valid_smiles(smiles):
            return jsonify({'success': False, 'error': 'Invalid chemical structure'})
        
        # Calculate properties
        properties = property_calculator.calculate_properties(smiles)
        
        return jsonify({
            'success': True,
            'smiles': smiles,
            'original_input': user_input,
            'input_type': input_type,
            'properties': properties
        })
        
    except Exception as e:
        app.logger.error(f"Error in calculate_properties: {str(e)}")
        return jsonify({'success': False, 'error': str(e)})

@app.route('/api/random_molecule')
def random_molecule():
    """Generate a random molecule for demonstration"""
    try:
        result = molecular_generator.generate_random()
        
        if result['success']:
            smiles = result['smiles']
            properties = property_calculator.calculate_properties(smiles)
            
            return jsonify({
                'success': True,
                'smiles': smiles,
                'properties': properties
            })
        else:
            return jsonify({
                'success': False,
                'error': result.get('error', 'Random generation failed')
            })
            
    except Exception as e:
        app.logger.error(f"Error in random_molecule: {str(e)}")
        return jsonify({'success': False, 'error': str(e)})

@app.route('/api/optimize_properties', methods=['POST'])
def optimize_properties():
    """Optimize molecule generation for target properties"""
    try:
        data = request.get_json() if request.is_json else request.form
        
        target_properties = {
            'qed': float(data.get('target_qed', 0.5)),
            'sa_score': float(data.get('target_sa', 3.0)),
            'logp': float(data.get('target_logp', 2.0))
        }
        
        # Run optimization (simplified approach)
        result = molecular_generator.optimize_for_properties(target_properties)
        
        if result['success']:
            molecules = []
            for smiles in result['molecules']:
                properties = property_calculator.calculate_properties(smiles)
                molecules.append({
                    'smiles': smiles,
                    'properties': properties
                })
            
            return jsonify({
                'success': True,
                'molecules': molecules,
                'optimization_info': result.get('info', {})
            })
        else:
            return jsonify({
                'success': False,
                'error': result.get('error', 'Optimization failed')
            })
            
    except Exception as e:
        app.logger.error(f"Error in optimize_properties: {str(e)}")
        return jsonify({'success': False, 'error': str(e)})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
