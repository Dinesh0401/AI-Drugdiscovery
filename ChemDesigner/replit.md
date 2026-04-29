# Molecular Discovery Platform

## Overview

This is an AI-powered molecular generation and analysis platform that combines graph neural networks with reinforcement learning for drug discovery applications. The system allows users to generate novel molecular structures by exploring latent space parameters and analyze existing molecules for their drug-like properties. The platform features an interactive web interface built with Flask and uses RDKit for molecular property calculations.

## User Preferences

Preferred communication style: Simple, everyday language.

## System Architecture

### Frontend Architecture
- **Template Engine**: Jinja2 templates with Bootstrap 5 for responsive UI
- **Styling**: Custom CSS with dark theme optimized for scientific applications
- **JavaScript**: Vanilla JavaScript modules for molecular visualization and user interactions
- **Canvas Rendering**: HTML5 Canvas-based molecular structure renderer with 2D visualization capabilities
- **Component Structure**: Modular template inheritance with base layout and specialized pages for generation and analysis

### Backend Architecture
- **Web Framework**: Flask application with RESTful API endpoints
- **Service Layer**: Modular service architecture with separate concerns:
  - `MolecularService`: SMILES validation and basic molecular operations
  - `PropertyCalculator`: Comprehensive molecular property calculations
  - `MolecularGenerator`: Simplified molecular generation from latent space parameters
- **Error Handling**: Comprehensive logging and graceful fallbacks when RDKit is unavailable
- **Session Management**: Flask sessions with configurable secret key

### Data Storage Solutions
- **In-Memory Storage**: Current implementation uses in-memory data structures for generated molecules
- **File-Based Assets**: Static molecular fragments and drug-like molecule templates stored in service classes
- **No Database**: Simplified architecture without persistent storage for rapid prototyping

### Authentication and Authorization
- **No Authentication**: Current implementation is open-access for demonstration purposes
- **Session Security**: Basic Flask session management with configurable secret key

## External Dependencies

### Core Scientific Libraries
- **RDKit**: Primary molecular informatics toolkit for SMILES validation, property calculations, and molecular descriptors
- **PyTorch Geometric**: Referenced for graph-based molecular representation (future implementation)
- **Fallback Systems**: Custom implementations when RDKit is unavailable to ensure platform functionality

### Web Framework Dependencies
- **Flask**: Core web framework for backend API and template rendering
- **Bootstrap 5**: Frontend CSS framework with dark theme customization
- **Font Awesome**: Icon library for user interface elements

### Property Calculation Features
- **Molecular Descriptors**: Molecular weight, LogP, hydrogen bond donors/acceptors
- **Drug-Likeness Metrics**: Lipinski's Rule of Five compliance
- **Synthesis Accessibility**: Synthetic accessibility scoring
- **QED Score**: Quantitative estimate of drug-likeness

### Development and Deployment
- **Logging**: Python logging module for debugging and monitoring
- **Environment Variables**: Configuration through environment variables for session secrets
- **Static Assets**: CSS and JavaScript served through Flask's static file handling