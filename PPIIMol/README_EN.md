# PPIIMoL – Automated Detection of PPII Helices in PyMOL

**PPIIMoL** is a Python module designed to integrate directly into PyMOL for the automatic detection of Polyproline II (PPII) helices in protein structures.

PPII helices, often present in glycine- and proline-rich proteins and relevant in neurobiological processes, are typically not explicitly annotated in PDB files. This tool identifies geometric patterns compatible with PPII helices and displays them directly in PyMOL, streamlining structural analysis.

## 🔬 Key Features

- Load local PDB files or fetch structures directly within PyMOL.
- Automatic removal of solvents and hydrogen addition.
- Calculation of φ and ψ angles for each residue.
- Detection of segments consistent with PPII conformation.
- Direct visualization of results in PyMOL using pseudoatoms.
- Export of detected segments to `.csv` and `.pdb` formats for external analysis.

## 📦 Requirements

- [PyMOL](https://pymol.org/) (version with Python scripting support)
- Python 3.8+
- Standard Python libraries (`math`, `tkinter`, `csv`, `os`)

## 🚀 Usage

1. Launch PyMOL with Python scripting support.
2. Run the script `PPIIMoL.py` via the PyMOL command line or graphical interface.
3. Load or fetch a protein structure.
4. Add hydrogens and perform the analysis.
5. Visualize the detected PPII helices and/or export the results.

## 🧪 Test Case

The module has been validated using the 3BOG protein, a well-known model with a high presence of PPII helices, confirming its accuracy and utility for structural research.

## 📄 License

Released under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).

## ✍️ Author

This module was developed as part of a Bachelor's Thesis in Computer Engineering (UNIR), in collaboration with the Protein Structure, Dynamics and Interactions Group at the Institute of Physical Chemistry “Blas Cabrera” (IQF-CSIC).

Author: Silvia Enma Rodríguez Fernández  
GitHub: [@silviaenma](https://github.com/silviaenma)

---

*Suggestions, contributions, and improvements are welcome.*
