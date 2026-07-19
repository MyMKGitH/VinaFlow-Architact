# 🧬 VinaFlow Architect (v2.1.0-Stable)

### **Automated High-Throughput Molecular Docking & Structural Analysis Pipeline**

VinaFlow Architect is an open-source, privacy-focused, cross-platform web application designed for automated structural biology workflows, molecular docking, and interactive virtual screening. 

By decoupling the computational core from rigid, platform-dependent source compilation, **Version 2.1.0** introduces a robust, native subprocess CLI routing layer. This architecture allows seamless high-throughput simulations natively across both **Windows and Linux environments** without requiring localized C++ compiler chains (e.g., Boost library dependencies).

---

## 🚀 Core Features

* **Multi-Mode Ligand Architectural Engine:** Automatically fetches, cleans, and prepares small molecules via PubChem API integration, processes custom SMILES strings, or handles automated peptide parsing.
* **Intelligent Search Optimization Strategies:**
    * *Global Blind Docking:* Automated whole-protein bounding box calculation.
    * *Heuristic Binding Pocket Prediction:* Active-site targeting based on protein geometric density.
    * *Manual Targeted Coordinates:* Precise grid control for refined docking parameters.
* **Dynamic RMSD Engagement Layer:** * Automatically executes heavy-atom RMSD calculation against reference structures during **Crystallographic Redocking Verification**.
    * Evaluates geometric structural drift and pose clustering against the highest-affinity model (Pose 1) during *de novo* screening campaigns.
* **Interactive 3D Complex Mapping:** Implements real-time 3D WebGL ribbon and ligand atom-level visualization directly within the dashboard using `stmol` and `py3Dmol`.
* **Privacy-First & Cloud-Ready:** Built to run completely offline on local machines or deploy seamlessly to containerized cloud architectures (such as Hugging Face Spaces).

---

## 🛠️ Tech Stack & Infrastructure

* **User Interface:** Streamlit (Python Dashboard Framework)
* **Cheminformatics & Geometry:** RDKit, Meeko, Gemmi, SciPy
* **Computational Core:** AutoDock Vina (CLI Executable Integration)
* **Visualization:** WebGL, py3Dmol, stmol

---

## 📦 Installation & Local Setup

### 1. Clone the Repository
```bash
git clone [https://github.com/MyMKGitH/VinaFlow-Architect.git](https://github.com/MyMKGitH/VinaFlow-Architect.git)
cd VinaFlow-Architect

---

[![DOI](https://zenodo.org/badge/1276759116.svg)](https://doi.org/10.5281/zenodo.21442287)
