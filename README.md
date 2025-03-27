# connexinmd.py
Connexin Molecular Dynamics Simulation Pipeline
End-to-End Pipeline for Modeling Cx34.7 and Cx35 Connexin Interactions
I've created a comprehensive computational pipeline to model and analyze wildtype Cx34.7 and Cx35 connexin interactions.
Overview of the Pipeline
The pipeline consists of the following components:

Main Pipeline Script (connexin-md-pipeline): Handles the end-to-end workflow from sequence acquisition through simulation to energy analysis.
NAMD Configuration (namd-config): Contains the molecular dynamics parameters for running simulations with NAMD.
Energy Analysis Script (energy-analysis): A VMD/NAMD-Energy script for calculating non-bonded interaction energies at the docking interface.
Visualization Script (visualization-script): Creates publication-quality visualizations of the energy analysis results.
Run Script (run-pipeline): A bash script to execute the entire pipeline from start to finish.

Key Features

Homology Modeling: Uses SWISS-MODEL to create structures of Cx34.7 and Cx35 based on available cryo-EM templates
System Setup: Embeds connexons in POPC membranes with appropriate ion concentrations
Simulation Protocol: Includes minimization, equilibration, and production phases
Energy Analysis: Calculates total docking energies and per-residue contributions
Visualization: Generates comprehensive visualizations including energy bar charts, residue energy profiles, and interaction heatmaps

How to Use

Setup: Install required dependencies (Python libraries, NAMD, VMD)
Run: Execute the run-pipeline.sh script
Analysis: The results will be in the connexin_md/analysis directory

Scientific Approach
The code implements the same approach as in your document:

Obtaining protein sequences from UniProtKB
Creating homology models based on high-resolution cryo-EM structures
Placing the models in POPC membranes with appropriate ion concentrations
Running minimization (2 ns) followed by equilibration (20 ns) and production (30 ns)
Calculating non-bonded interaction energies at the docking interface

Expected Results
The pipeline will generate comprehensive data on the interaction energies of:

Homotypic Cx34.7-Cx34.7 channels
Homotypic Cx35-Cx35 channels
Heterotypic Cx34.7-Cx35 channels

Based on similar analyses in the literature, we would expect homotypic channels to have more favorable (more negative) interaction energies than heterotypic channels, providing insights into connexin compatibility and gap junction formation preferences.
