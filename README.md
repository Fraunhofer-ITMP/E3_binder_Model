# <center> Pharmacophore-based ML model to predict ligand selectivity for E3 ligase binders </center>

Model based on E3Binder list from Protac-DB/ProtacPedia and Evolvus. 

Under `data` folder the three data files containing structural details of molecules used for selectivity modeling.

* **Asinex library** can be found also [here](https://www.asinex.com/protein-degradation).
* **Broad library version 2018** used can be found [here](https://clue.io/repurposing#download-data).

**R script** has been setup with RStudio v.1.21335 and R v.3.6.1 All libraries used in R are cited in the script. This script can be found under the `script` directory.

Due to licensing restrictions of Evolvus data, we have removed smiles notation from the file. However, ErG description can be still used to rerun script and predictions.
**ErG Script** can be found in original format as MOE SVL script [here](https://svl.chemcomp.com/data/Extended_Reduced_Graph__ErG__fingerprint.svlx) and searching for keyword "Reduced" or as RDkit fingeprint scheme [here](https://www.rdkit.org/docs/GettingStartedInPython.html)

> NOTE: Manuscript submitted. Pre-print [available](https://chemrxiv.org/engage/chemrxiv/article-details/643d416808c86922ff297ab9)!
