This code contains a collection of python scripts for calculating structural properties of protein complexes.
 
# Install
Download the code and run
```
python setup.py install
```

# How to execute the script
1. Create a text file with the proteins of interest and the interacting chains where the ":" symbol indicates the protein-protein interface that is being investigated. For example:
  ```
  1AKJ_AB:DE
  1AK4_A:D
  ```
2. Invoke setupPpiDb.py on the file you just created. For example (assuming it is saved as PDBs.txt):
  ```
  python setupPpiDb.py PDBs.txt
  ```
Follow the script instructions and that's it!

# PPI features
The above script will:
* Download PDBs from the RCSB PDB
* Add hydrogens using [molprobity](http://molprobity.biochem.duke.edu/)
* Calculate/identify the following properties:
  * Accessible surface area ([ASA](https://en.wikipedia.org/wiki/Accessible_surface_area))
  * interface atoms/residues in the molecule (ΔASA > 0 or Distance > 4)
  * periphery index
  * hydrogen bonds
  * Van der Waals
  * Electrostatic charges
* All the results are saved to csv files and to SQL database


# Credits
All the scripts were developed in [Dr. Julia Shifman lab](http://bio.huji.ac.il/shifman/index.html).

1. Erijman, Ariel, Eran Rosenthal, and Julia Shifman, [How Structure Defines Affinity in Protein-Protein Interactions](http://dx.doi.org/10.1371/journal.pone.0110085). PLOS one 9.10 (2014)
