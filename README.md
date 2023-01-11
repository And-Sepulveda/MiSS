# Micellar Supramolecular Structure Software (MiSS)*

Miss** is a software dedicated to the identification of supramolecular structures and calculation of lattice parameters, obtained from Small-Angle Scattering Neutrons (SANS) and X-ray (SAXS). The error is calculated from a Gaussian fit over the principal peak.

# Structures identified 
- Cubic Fm3m;
- Cubic Pm3n;
- Cubic Pn3m;
- Cubic Ia3d;
- Cubic Im3m;
- Centred face cubic;
- Centred body cubic;
- Lamellar and
- Hexagonal 

# How to use
- Download MiSS_2.1.py file;
- Open the file;
- The SANS or SAXS data should be saved in .txt extension;
- In MiSS.py, change the variable "arq" to file path where the data is saved;
- Variable "var" defines the distance between principal peak and parameter calculated
- Variable "div" defines how many divisions in x axis;
- To create Gaussian fit over the principal peak :
  - choose the first min peak (variable "curv_ini");
  - choose the second min peak (variable "curv_fin");
- Run the program.

# Notes:
- Software registered under Brazilian agency INPI (BR 51 2021 000746-8)
- If you use this software please cite "Sepulveda et al. Supramolecular structure organization and rheological properties modulate the performance of hyaluronic acid-loaded thermosensitive hydrogels as drug-delivery systems. Journal of Colloid and Interface Science 630 (2023) 328–340" 
