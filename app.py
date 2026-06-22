import streamlit as st
import os
import sys
import requests
import subprocess
import shutil
import platform
import numpy as np
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
import py3Dmol
from stmol import showmol

# Initialize application-wide file isolation
WORKING_DIR = "vinaflow_runtime"
os.makedirs(WORKING_DIR, exist_ok=True)
if not os.getcwd().endswith(WORKING_DIR):
    os.chdir(WORKING_DIR)

st.set_page_config(page_title="VinaFlow Architect v2.0.0", layout="wide")
st.title("🧬 VinaFlow Architect `v2.0.0` (Cross-Platform CLI Engine)")
st.caption("Automated High-Throughput Molecular Docking & Structural Analysis Pipeline")

# --- CROSS-PLATFORM BINARY ROUTER ---

def locate_vina_executable():
    """Dynamically routes paths to find local Windows binaries or Linux system calls."""
    # Step up a level if execution trapped inside runtime workspace
    parent_dir = os.path.dirname(os.getcwd())
    
    if platform.system() == "Windows":
        # Look for the file you downloaded locally
        for binary in ["vina.exe", "vina_1.2.7_win.exe"]:
            if os.path.exists(os.path.join(parent_dir, binary)):
                return os.path.abspath(os.path.join(parent_dir, binary))
            if os.path.exists(binary):
                return os.path.abspath(binary)
        return "vina.exe"  # System fallback assumed inside PATH environmental routes
    else:
        # Cloud Target Environment Routing (Hugging Face Linux Engines)
        for binary in ["autodock-vina", "vina"]:
            if shutil.which(binary):
                return binary
        return "vina"

# --- CORE COMPUTATIONAL ENGINE FUNCTIONS ---

def fetch_pdb_protein(pdb_id):
    """Downloads a raw PDB file directly from the RCSB Server."""
    url = f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        pdb_path = f"{pdb_id.lower()}.pdb"
        with open(pdb_path, "w") as f:
            f.write(response.text)
        return pdb_path
    return None

def prepare_protein_obabel(input_pdb, output_pdbqt):
    """Uses OpenBabel CLI to clean receptors, strip water/heteroatoms, and fix charges."""
    cmd = f"obabel {input_pdb} -O {output_pdbqt} -xr -d -h"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return os.path.exists(output_pdbqt)

def smiles_to_pdbqt(smiles, output_pdbqt):
    """Generates an optimized, energy-minimized 3D conformer from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    
    preparer = MoleculePreparation()
    preparer.prepare(mol)
    preparer.write_pdbqt_file(output_pdbqt)
    return os.path.exists(output_pdbqt)

def peptide_to_pdbqt(sequence, output_pdbqt):
    """Translates an amino acid sequence string into an energy-minimized 3D structural model."""
    mol = Chem.MolFromSequence(sequence)
    if not mol:
        return False
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    
    preparer = MoleculePreparation()
    preparer.prepare(mol)
    preparer.write_pdbqt_file(output_pdbqt)
    return os.path.exists(output_pdbqt)

def calculate_heavy_atom_rmsd(ref_pdbqt, target_pdbqt, pose_index=1):
    """Computes coordinate-based RMSD for validation without topological graph dependencies."""
    def parse_heavy_atoms(file_path, target_pose=1):
        coords = []
        current_pose = 1
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith("MODEL"):
                    current_pose = int(line.split()[1])
                if current_pose != target_pose:
                    continue
                if line.startswith(("ATOM", "HETATM")):
                    atom_type = line[76:78].strip()
                    if atom_type != "H" and atom_type != "":
                        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        return np.array(coords)

    ref_coords = parse_heavy_atoms(ref_pdbqt, target_pose=1)
    target_coords = parse_heavy_atoms(target_pdbqt, target_pose=pose_index)
    
    if len(ref_coords) == 0 or len(ref_coords) != len(target_coords):
        return None
    
    diff = ref_coords - target_coords
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

def run_heuristic_pocket_detection(pdb_path):
    """Heuristic Engine: Identifies coordinate density hotspots of aromatic & functional residues."""
    coords = []
    target_residues = ["HIS", "ASP", "SER", "CYS", "LYS", "TYR", "TRP", "PHE"]
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                res_name = line[17:20].strip()
                if res_name in target_residues:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not coords:
        return [0.0, 0.0, 0.0], [30.0, 30.0, 30.0]
    
    coords = np.array(coords)
    centroid = coords.mean(axis=0)
    std_dev = np.std(coords, axis=0)
    box_size = np.clip(std_dev * 2.5, 18.0, 32.0)
    return centroid.tolist(), box_size.tolist()

def get_blind_box_params(pdb_path):
    """Calculates geometric boundaries to cleanly encompass the global protein structure."""
    coords = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    coords = np.array(coords)
    min_c, max_c = coords.min(axis=0), coords.max(axis=0)
    center = ((min_c + max_c) / 2).tolist()
    size = ((max_c - min_c) + 12.0).tolist()
    return center, size

def convert_pdbqt_to_pdb(input_pdbqt, output_pdb):
    """Converts structural files back to clean PDB models for WebGL visualization engines."""
    cmd = f"obabel {input_pdbqt} -O {output_pdb}"
    subprocess.run(cmd, shell=True, capture_output=True)

def parse_vina_output_energies(output_pdbqt):
    """Parses structural PDBQT result text matrices to extract raw affinity values."""
    energies = []
    if not os.path.exists(output_pdbqt):
        return energies
    with open(output_pdbqt, "r") as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 4:
                    energies.append(float(parts[3]))
    return energies

# --- USER INTERFACE CONFIGURATION ---

st.sidebar.header("🛠️ Pipeline Parameters")
protein_mode = st.sidebar.selectbox("Target Receptor Structural Input", ["RCSB PDB Database Lookup"])
pdb_id = st.sidebar.text_input("Enter 4-Character PDB ID", value="1IEP")

ligand_mode = st.sidebar.selectbox("Ligand Architecture Input", [
    "PubChem Chemical Database Engine", 
    "De Novo Peptide Sequence Parser", 
    "Native Crystallographic Ligand (Redocking Verification)"
])

ligand_input = ""
if ligand_mode == "PubChem Chemical Database Engine":
    ligand_input = st.sidebar.text_input("Enter Compound Common Name or Canonical SMILES", value="Imatinib")
elif ligand_mode == "De Novo Peptide Sequence Parser":
    ligand_input = st.sidebar.text_input("Enter FASTA Single-Letter Peptide String", value="GAV")
elif ligand_mode == "Native Crystallographic Ligand (Redocking Verification)":
    ligand_input = st.sidebar.text_input("Enter 3-Letter Heteroatom Residue ID (e.g., STI, HEM)", value="STI")

st.sidebar.markdown("---")
docking_strategy = st.sidebar.radio("Grid Search Optimization Strategy", [
    "Global Blind Docking", 
    "Automated Binding Pocket Prediction (Heuristic Engine)",
    "Manual Targeted Coordinates (Active Site Coordinates)"
])

center_x, center_y, center_z = 0.0, 0.0, 0.0
size_x, size_y, size_z = 20.0, 20.0, 20.0

if docking_strategy == "Manual Targeted Coordinates (Active Site Coordinates)":
    st.sidebar.markdown("### Grid Center Target Matrix")
    center_x = st.sidebar.number_input("Center Vector X", value=0.0)
    center_y = st.sidebar.number_input("Center Vector Y", value=0.0)
    center_z = st.sidebar.number_input("Center Vector Z", value=0.0)
    st.sidebar.markdown("### Inner Search Volume Box Dimensions (Å)")
    size_x = st.sidebar.number_input("Dimension X (Å)", value=20.0)
    size_y = st.sidebar.number_input("Dimension Y (Å)", value=20.0)
    size_z = st.sidebar.number_input("Dimension Z (Å)", value=20.0)

exhaustiveness = st.sidebar.slider("Vina Calculation Exhaustiveness Engine", min_value=1, max_value=16, value=8)

# --- EXECUTION ENGINE FLOW ---

if st.sidebar.button("🚀 Execute Structural Docking"):
    if not pdb_id or not ligand_input:
        st.error("Halting: Structural target fields or ligand identities must be fully configured.")
        st.stop()

    vina_path = locate_vina_executable()

    # Step 1: Parse and Clean Receptors
    with st.spinner("Executing Pipeline Tier 1: Downloading & formatting target protein matrix..."):
        pdb_file = fetch_pdb_protein(pdb_id)
        if not pdb_file:
            st.error(f"Network Fault: Unable to acquire PDB structure matching ID: {pdb_id}")
            st.stop()
        
        receptor_pdbqt = "receptor.pdbqt"
        if not prepare_protein_obabel(pdb_file, receptor_pdbqt):
            st.error("System Error: OpenBabel structural normalization run failed.")
            st.stop()

    # Step 2: Compute Ligand Structures
    with st.spinner("Executing Pipeline Tier 2: Resolving, minimizing, and converting ligand geometries..."):
        ligand_pdbqt = "ligand.pdbqt"
        native_validated_pdbqt = "native_ligand.pdbqt"
        success = False
        
        if ligand_mode == "PubChem Chemical Database Engine":
            if "=" in ligand_input or (len(ligand_input) > 12 and any(char in ligand_input for char in ['@', '(', ')'])):
                smiles = ligand_input
            else:
                try:
                    compounds = pcp.get_compounds(ligand_input, 'name')
                    if compounds:
                        smiles = compounds[0].canonical_smiles
                    else:
                        st.error("Database Lookup Failure: Compound identity not recognized inside PubChem repositories.")
                        st.stop()
                except Exception as error:
                    st.error(f"PubChem Query Failure: {error}")
                    st.stop()
            success = smiles_to_pdbqt(smiles, ligand_pdbqt)
            
        elif ligand_mode == "De Novo Peptide Sequence Parser":
            success = peptide_to_pdbqt(ligand_input, ligand_pdbqt)
            
        elif ligand_mode == "Native Crystallographic Ligand (Redocking Verification)":
            native_raw_pdb = "native_extracted.pdb"
            with open(pdb_file, "r") as incoming_pdb, open(native_raw_pdb, "w") as target_out:
                for line in incoming_pdb:
                    if line.startswith("HETATM") and ligand_input.upper() in line:
                        target_out.write(line)
            
            if os.path.getsize(native_raw_pdb) == 0:
                st.error(f"Structural Context Error: Structural code '{ligand_input.upper()}' missing inside targeted PDB record.")
                st.stop()
                
            cmd = f"obabel {native_raw_pdb} -O {ligand_pdbqt} -h"
            subprocess.run(cmd, shell=True)
            subprocess.run(f"obabel {native_raw_pdb} -O {native_validated_pdbqt} -h", shell=True)
            success = os.path.exists(ligand_pdbqt)

        if not success:
            st.error("Execution Interrupted: Failure encountered during ligand structural parsing operations.")
            st.stop()

    # Step 3: Space Matrix Computations
    if docking_strategy == "Global Blind Docking":
        center, size = get_blind_box_params(pdb_file)
    elif docking_strategy == "Automated Binding Pocket Prediction (Heuristic Engine)":
        center, size = run_heuristic_pocket_detection(pdb_file)
    else:
        center = [center_x, center_y, center_z]
        size = [size_x, size_y, size_z]

    # Step 4: Run AutoDock Vina CLI via Subprocess
    with st.spinner("Executing Pipeline Tier 3: Invoking AutoDock Vina Engine Simulation..."):
        output_poses_pdbqt = "docked_output_poses.pdbqt"
        
        # Build argument parameters array for standard CLI execution
        vina_cmd = [
            vina_path,
            "--receptor", receptor_pdbqt,
            "--ligand", ligand_pdbqt,
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
            "--exhaustiveness", str(exhaustiveness),
            "--out", output_poses_pdbqt
        ]
        
        try:
            # Trigger OS execution profile safely
            run_result = subprocess.run(vina_cmd, capture_output=True, text=True, check=True)
            
            # Parse output matrix
            energy_scores = parse_vina_output_energies(output_pdbqt=output_poses_pdbqt)
            
            if not energy_scores:
                st.error("Simulation Warning: Structural poses were not generated. Verify site coordinate configurations.")
                st.stop()
                
            col1, col2 = st.columns([3, 2])
            
            with col1:
                st.subheader("📊 Docking Run Final Score Matrices")
                results_table = []
                for idx, delta_g in enumerate(energy_scores):
                    kd_value = np.exp(delta_g / 0.5925)
                    if kd_value < 1e-6:
                        kd_display = f"{kd_value * 1e9:.2f} nM"
                    elif kd_value < 1e-3:
                        kd_display = f"{kd_value * 1e6:.2f} µM"
                    else:
                        kd_display = f"{kd_value * 1e3:.2f} mM"
                    
                    row = {
                        "Pose Index": idx + 1,
                        "Estimated ΔG (kcal/mol)": round(delta_g, 2),
                        "Affinity Constant (Est. Kd)": kd_display
                    }
                    
                    if ligand_mode == "Native Crystallographic Ligand (Redocking Verification)":
                        calculated_rmsd = calculate_heavy_atom_rmsd(native_validated_pdbqt, output_poses_pdbqt, pose_index=idx+1)
                        row["Validation Heavy-Atom RMSD (Å)"] = round(calculated_rmsd, 2) if calculated_rmsd is not None else "Error"
                        
                    results_table.append(row)
                
                st.table(results_table)
                
                st.subheader("💾 Export Calculations")
                with open(output_poses_pdbqt, "rb") as final_f:
                    st.download_button(
                        label="📥 Download Output Coordinates (PDBQT Format)",
                        data=final_f,
                        file_name=f"VinaFlow_{pdb_id}_poses.pdbqt",
                        mime="text/plain"
                    )

            with col2:
                st.subheader("🔮 Interactive 3D Complex Map")
                with st.spinner("Rendering WebGL Model Viewer..."):
                    receptor_view_pdb = "receptor_display.pdb"
                    poses_view_pdb = "poses_display.pdb"
                    convert_pdbqt_to_pdb(receptor_pdbqt, receptor_view_pdb)
                    convert_pdbqt_to_pdb(output_poses_pdbqt, poses_view_pdb)
                    
                    if os.path.exists(receptor_view_pdb) and os.path.exists(poses_view_pdb):
                        webgl_viewer = py3Dmol.view(width=550, height=500)
                        webgl_viewer.addModel(open(receptor_view_pdb, 'r').read(), 'pdb')
                        webgl_viewer.setStyle({'cartoon': {'color': 'spectrum'}})
                        
                        webgl_viewer.addModel(open(poses_view_pdb, 'r').read(), 'pdb')
                        webgl_viewer.setStyle({'model': 1}, {'stick': {'colorscheme': 'cyanCarbon'}})
                        
                        webgl_viewer.zoomTo()
                        showmol(webgl_viewer, height=500, width=550)
                    else:
                        st.info("System Notification: Local formatting limits prevented 3D execution.")
                        
        except subprocess.CalledProcessError as e:
            st.error("System Execution Blocked: Vina binary failed to compute spatial properties.")
            st.text(e.stderr)
        except Exception as runtime_error:
            st.error(f"Structural Simulation Failure: {runtime_error}")
