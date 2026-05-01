import json
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _build_wsl_python_script(bond_lengths: list[float]) -> str:
    grid_json = json.dumps(bond_lengths)
    return f"""
import json

from pyscf import fci, gto, scf
from qiskit.primitives import StatevectorEstimator
from qiskit_algorithms.minimum_eigensolvers import VQE
from qiskit_algorithms.optimizers import SLSQP
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.units import DistanceUnit


def pyscf_reference_energy(r_angstrom: float) -> float:
    atom = f"H 0 0 0; H 0 0 {{r_angstrom}}"
    mol = gto.M(atom=atom, basis="sto3g", unit="Angstrom", spin=0, charge=0, verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    cisolver = fci.FCI(mf)
    e_fci, _ = cisolver.kernel()
    return float(e_fci)


def qiskit_vqe_energy(r_angstrom: float) -> float:
    atom = f"H 0 0 0; H 0 0 {{r_angstrom}}"
    driver = PySCFDriver(atom=atom, basis="sto3g", charge=0, spin=0, unit=DistanceUnit.ANGSTROM)
    problem = driver.run()

    mapper = JordanWignerMapper()
    qubit_op = mapper.map(problem.hamiltonian.second_q_op())

    init_state = HartreeFock(problem.num_spatial_orbitals, problem.num_particles, mapper)
    ansatz = UCCSD(problem.num_spatial_orbitals, problem.num_particles, mapper, initial_state=init_state)
    vqe = VQE(StatevectorEstimator(), ansatz, SLSQP(maxiter=120))

    result = vqe.compute_minimum_eigenvalue(qubit_op)
    interpreted = problem.interpret(result)
    return float(interpreted.total_energies[0].real)


bond_lengths = {grid_json}
reference_curve = []
vqe_curve = []

for r in bond_lengths:
    e_ref = pyscf_reference_energy(r)
    e_vqe = qiskit_vqe_energy(r)
    reference_curve.append(e_ref)
    vqe_curve.append(e_vqe)

print(json.dumps({{"bond_lengths": bond_lengths, "reference": reference_curve, "vqe": vqe_curve}}))
""".strip()


def run_wsl_h2_curves(bond_lengths: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    script_body = _build_wsl_python_script([float(x) for x in bond_lengths])
    temp_script = Path("_wsl_h2_curve_runner.py")
    temp_script.write_text(script_body, encoding="utf-8")

    command = [
        "wsl",
        "-e",
        "bash",
        "-lc",
        "cd '/mnt/e/QC Project/Implementation/Molecular Binding using VQE' "
        "&& . .wsl_venv/bin/activate "
        "&& python _wsl_h2_curve_runner.py",
    ]

    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        if completed.returncode != 0:
            raise RuntimeError(
                "WSL calculation failed.\n"
                f"STDOUT:\n{completed.stdout}\n\nSTDERR:\n{completed.stderr}"
            )

        lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
        payload = json.loads(lines[-1])
        return (
            np.array(payload["bond_lengths"], dtype=float),
            np.array(payload["reference"], dtype=float),
            np.array(payload["vqe"], dtype=float),
        )
    finally:
        if temp_script.exists():
            temp_script.unlink()


def main() -> None:
    print("Running H2 binding-curve comparison (PySCF reference vs Qiskit VQE simulator) ...")
    bond_grid = np.array(sorted(set(np.round(np.linspace(0.4, 1.6, 13), 2).tolist() + [0.74])), dtype=float)
    bond_lengths, reference_energies, vqe_energies = run_wsl_h2_curves(bond_grid)

    abs_error = np.abs(reference_energies - vqe_energies)
    ref_min_idx = int(np.argmin(reference_energies))
    vqe_min_idx = int(np.argmin(vqe_energies))

    print("\nCurve summary:")
    print(f"Reference minimum: R = {bond_lengths[ref_min_idx]:.2f} A, E = {reference_energies[ref_min_idx]:.8f} Ha")
    print(f"VQE minimum      : R = {bond_lengths[vqe_min_idx]:.2f} A, E = {vqe_energies[vqe_min_idx]:.8f} Ha")
    print(f"Max absolute error: {abs_error.max():.3e} Ha")
    print(f"Mean absolute error: {abs_error.mean():.3e} Ha")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9), sharex=True)

    ax1.plot(bond_lengths, reference_energies, marker="o", linewidth=2.0, label="Original/reference (PySCF FCI)")
    ax1.plot(bond_lengths, vqe_energies, marker="s", linestyle="--", linewidth=2.0, label="Simulator (Qiskit VQE)")
    ax1.axvline(0.74, color="#888888", linestyle=":", linewidth=1.5, label="Experimental reference ~0.74 A")
    ax1.set_ylabel("Total Energy (Hartree)")
    ax1.set_title("H2 Potential Energy vs Bond Length")
    ax1.grid(alpha=0.25)
    ax1.legend()

    ax2.plot(bond_lengths, abs_error, marker="d", color="#d62728", linewidth=1.8)
    ax2.set_xlabel("Bond Length (A)")
    ax2.set_ylabel("Absolute Error (Hartree)")
    ax2.set_title("|Reference - VQE| Across Bond Lengths")
    ax2.grid(alpha=0.25)

    plt.tight_layout()

    output_path = Path("h2_pe_curve_comparison.png")
    plt.savefig(output_path, dpi=220)

    if "agg" in plt.get_backend().lower():
        print(f"\nSaved plot to {output_path}")
    else:
        plt.show()


if __name__ == "__main__":
    main()