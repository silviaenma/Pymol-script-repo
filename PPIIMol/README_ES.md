# PPIIMoL - Detección automática de hélices PPII en PyMOL

**PPIIMoL** es un módulo en Python diseñado para integrarse directamente con PyMOL y automatizar la detección de hélices de Poliprolina II (PPII) en estructuras proteicas.

Estas hélices, relevantes en procesos neurobiológicos y en proteínas ricas en glicina y prolina, no suelen estar identificadas explícitamente en los archivos PDB. Este módulo permite detectar automáticamente patrones geométricos compatibles con hélices PPII y visualizarlos en PyMOL, facilitando su análisis estructural.

## 🔬 Funcionalidades principales

- Carga de archivos PDB o descarga directa desde PyMOL.
- Eliminación automática de solventes y adición de hidrógenos.
- Cálculo de ángulos φ y ψ para cada residuo.
- Detección de segmentos con conformación PPII.
- Visualización directa en PyMOL mediante pseudóatomos.
- Exportación de resultados a `.csv` y `.pdb` para análisis externo.

## 📦 Requisitos

- [PyMOL](https://pymol.org/) (versión con soporte de scripts Python)
- Python 3.8+
- Módulos estándar de Python (`math`, `tkinter`, `csv`, `os`)

## 🚀 Cómo usar

1. Abre PyMOL con soporte de scripts.
2. Ejecuta el módulo `PPIIMoL.py` desde la consola de PyMOL o usa el GUI incluido.
3. Carga o descarga una proteína.
4. Añade hidrógenos y ejecuta el análisis.
5. Visualiza las hélices PPII detectadas y/o exporta los datos.

## 🧪 Caso de prueba

La herramienta ha sido validada utilizando la proteína 3BOG, una estructura rica en hélices PPII, demostrando su eficacia y utilidad para análisis sistemático.

## 📄 Licencia

Publicado bajo licencia [GPL v3.0](https://www.gnu.org/licenses/gpl-3.0.html).

## ✍️ Autoría

Este módulo ha sido desarrollado como parte de un Trabajo de Fin de Grado en Ingeniería Informática (UNIR), en colaboración con el grupo de investigación del Instituto de Química-Física “Blas Cabrera” (IQF-CSIC).

Autora: Silvia Enma Rodríguez Fernández  
GitHub: [@silviaenma](https://github.com/silviaenma)

---

*Contribuciones, sugerencias y mejoras son bienvenidas.*
