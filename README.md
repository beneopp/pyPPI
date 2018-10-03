This is collection of simple python scripts for calculating basic structural properties/features of protein complexes.
 
# Install
Download the code and run
```
python setup.py install
```

# Usage
You can either use the scripts directly or use master script setupPpiDb to run most of them on a list of proteins.

For basic usage:
1. Create a text file with list of proteins and the interacting chains, for example:
  ```
  1AKJ_AB:DE
  1AK4_A:D
  ```
2. invoke setupPpiDb.py on the file you just created, for example (assuming it is saved as PDBs.txt) 
   ```
  python setupPpiDb.py PDBs.txt
  ```
Follow the script instructions and that's it!

# PPI features
The above script will:
* Download PDBs from the RCSB PDB
* Add hydrogens using [molprobity](http://molprobity.biochem.duke.edu/)
* Calculate various properties:
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