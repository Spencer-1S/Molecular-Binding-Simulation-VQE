# Hydrogen Molecule Binding Simulation using VQE

This project implements an H2 potential-energy curve workflow with:

- Original/reference curve from first-principles electronic structure in PySCF (run in WSL).
- Simulator curve from Qiskit VQE on the same molecular Hamiltonian.
- Side-by-side comparison plots: energy curves and absolute error vs bond length.

## Files

- [H2_VQE_Binding.ipynb](H2_VQE_Binding.ipynb): notebook walkthrough with explanations, curve comparison, and metrics.
- [Main.py](Main.py): script that runs WSL-based calculations and saves comparison plots.

## Setup

1) Windows environment (for plotting and orchestration):

```powershell
python -m venv .venv
& .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install numpy matplotlib
```

2) WSL environment (for PySCF + Qiskit chemistry calculations):

```powershell
wsl -e bash -lc "sudo apt-get update && sudo apt-get install -y python3.12-venv"
wsl -e bash -lc "cd '/mnt/e/QC Project/Implementation/Molecular Binding using VQE' && python3 -m venv .wsl_venv"
wsl -e bash -lc "cd '/mnt/e/QC Project/Implementation/Molecular Binding using VQE' && . .wsl_venv/bin/activate && python -m pip install --upgrade pip"
wsl -e bash -lc "cd '/mnt/e/QC Project/Implementation/Molecular Binding using VQE' && . .wsl_venv/bin/activate && python -m pip install pyscf qiskit qiskit-aer qiskit-algorithms qiskit-nature"
```

## Run the notebook

Open [H2_VQE_Binding.ipynb](H2_VQE_Binding.ipynb) in VS Code and run the cells from top to bottom.

The notebook includes:

- environment checks
- WSL execution helper
- reference vs simulator PE curve comparison
- absolute error analysis and minimum-energy location reporting

## Run the script

```powershell
& .venv\Scripts\python.exe Main.py
```

The script saves [h2_pe_curve_comparison.png](h2_pe_curve_comparison.png) in headless mode.

## Notes

- Energies are reported in Hartree.
- Bond length is reported in Angstrom.
- The experimental H2 equilibrium distance is approximately 0.74 Angstrom; your sampled-grid minimum should appear near this value.