# 🧬 PPIIMoL

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue)](https://www.python.org/)
[![PyMOL](https://img.shields.io/badge/PyMOL-2.x-green)](https://pymol.org/)

## 📖 Descripción

**PPIIMoL** es un módulo en Python diseñado para integrarse con **PyMOL** y automatizar la detección de hélices de poliprolina II (PPII) en proteínas. Esta herramienta facilita el análisis estructural, identificando rápidamente ángulos torsionales (`phi` y `psi`) característicos, así como posibles enlaces de hidrógeno no canónicos (Cα-H···O=C).

Desarrollado como parte de un **Trabajo de Fin de Grado en Ingeniería Informática** en colaboración con el laboratorio de neurociencia del CSIC.

---

## 🚀 Características
- 🔍 Detección automática de segmentos PPII mediante análisis de ángulos phi y psi.
- 🧬 Identificación de interacciones Cα-H···O=C relevantes para la estabilidad estructural.
- 📊 Exportación de resultados en CSV para análisis adicionales.
- 🎨 Visualización directa en PyMOL con códigos de color personalizados.
- 🖱️ Interfaz gráfica sencilla basada en Tkinter.

---

## 🛠️ Requisitos
- Python >= 3.9
- [PyMOL](https://pymol.org/2/) (entorno gráfico)
- Tkinter (incluido en la mayoría de distribuciones de Python)

---

## 📦 Instalación

1. Clona este repositorio:
   ```bash
   git clone https://github.com/silviaenma/PPIIMoL.git

Abre PyMOL y añade el módulo a la ruta de plugins o cárgalo manualmente:
run PPIIMoL/PPIIMoL.py

🧪 Ejemplo de uso
# Cargar el módulo en PyMOL
run PPIIMoL/PPIIMoL.py

# Detectar hélices PPII en un archivo PDB
load 3bog.pdb
ppii_detect()

📂 Los resultados se exportarán automáticamente en una carpeta organizada con fecha.

📜 Licencia
Este proyecto está bajo la licencia GNU GPLv3.

🌐 Description (English)
PPIIMoL is a Python module designed to integrate with PyMOL for the automatic detection of polyproline II (PPII) helices in proteins. This tool simplifies structural analysis by quickly identifying characteristic torsional angles (phi and psi), as well as potential non-canonical hydrogen bonds (Cα-H···O=C).

Developed as part of a Bachelor's Thesis in Computer Engineering in collaboration with the CSIC neuroscience lab.

🚀 Features
🔍 Automatic detection of PPII segments via phi/psi angle analysis.

🧬 Identification of Cα-H···O=C interactions relevant to structural stability.

📊 Export of results in CSV for further analysis.

🎨 Direct visualization in PyMOL with customizable color codes.

🖱️ Simple GUI based on Tkinter.

🛠️ Requirements
Python >= 3.9

PyMOL (graphical environment)

Tkinter (usually included in Python distributions)

📦 Installation
1. Clone this repository:
git clone https://github.com/silviaenma/PPIIMoL.git

2. Open PyMOL and load the module:
run PPIIMoL/PPIIMoL.py

🧪 Example usage
python
Copiar
Editar
# Load the module in PyMOL
run PPIIMoL/PPIIMoL.py

# Detect PPII helices in a PDB file
load 3bog.pdb
ppii_detect()
📂 Results are automatically exported in a dated folder.

📜 License
This project is licensed under the GNU GPLv3.

