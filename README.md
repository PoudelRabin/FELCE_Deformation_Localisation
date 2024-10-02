# FELCE_Deformation_Localisation
The repository contains finite element simulation code assocaited with following paper:

"Defomation Localisation in stretched liquid crystal elastomers"

Authors: Rabin Poudel, Yasemin Sengul, Angela Mihai

**Note**: This repository is for hosting code of the published paper and is not intended for regular changes and updates. The updated code if available will be available in github repository [FELCE_Deformation_Localisation_Updated](https://github.com/PoudelRabin/FELCE_Deformation_Localisation_Updated/)
## Folders
### src
This folder contains the sorce code for the finite element simualtion. The simulations are done in the open software [FEBio](https://febio.org/) ([FEBio Github](https://github.com/febiosoftware/FEBio)). FEBio 4.4 is used for the simualtion. Most of the codes are the FEBio 4.4 with user codes written by authors starts with the 'FELCE' and are in the folder src/FEBioMech. Some of the original files from src/FEBioMech are edited to register the FELCE codes.

### FEMaterialModel
This folder contains the two folders LCENecking and LCEStriping containing containing the material code for the necking and striping phenomena described in the paper. Each folders contain the c++ file "FELCEMaterial444" which need to be copied, moved to src/FEBioMech replacing the original according the phenomena that need to be simulated. 

### Examples
This folder contains the input file for simualting the LCE material. README file in the folder describes model.

## Building and running the model
Building and running is similar to the one building the original FEBio excutables. Please follow the README.md and BUILD.md in the src folder as well as read the documentation provided. More resources are availabe in [FEBio](https://febio.org/) ([FEBio Github](https://github.com/febiosoftware/FEBio)). Please note that software requres the pardiso solver for solving the equations. Please install intel oneAPI math kernel library.

