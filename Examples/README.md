## Model description

### model-1
This contains the input file for running the simulation described in fig 6 in the paper mentioned. This is the necking phenomena and the material model from the folder LCEMaterialModel/LCENecking need to be used.

### model-2 
This contains the input file for running the simulation described in fig 7 (last row with L=0.4) in the paper mentioned. This is the striping phenomena and the material model from the folder LCEMaterialModel/LCEStriping need to be used.

### model-3
This contains the input file for running the simulation described in fig 7 (first row with L=0.1) and fig 8 ((first row with L=0.1)  in the paper mentioned. This is the striping phenomena and the material model from the folder LCEMaterialModel/LCEStriping need to be used.

### model-4
This contains the input file for running the simulation described  fig 8 (second row with L=0.1, when frank constant increased by 100)  in the paper mentioned. This is the striping phenomena and the material model from the folder LCEMaterialModel/LCEStriping need to be used.

### model-5
This contains the input file for running the simulation described  fig 8 (second row with L=0.1, when frank constant increased by 1000)  in the paper mentioned. This is the striping phenomena and the material model from the folder LCEMaterialModel/LCEStriping need to be used.

## Material parameters
The variable used for the material parameters in input files and code.

**mu**: Found in the input file and resembles the $\mu_2$ used in the paper. It is the shear modulus of the neo-classical free energy. It can be changed in the input file.

**K**: Frank constant denoted by $K$ in the paper.

**a**: Stretch parameter denoted by $a$ in the paper. It can be updated in the model for striping phenomenon. Please note that for necking phenomena, $a$ is the particular function of the strain and is defined in the code itself. So, change in the value of $a$ in the input file for necking phenomena would not have any effects.

**B**: It is the bulk modulus needed for implementing nearly incompressible condition. Usually $B = 1000(\mu_1+\mu_2)$.

**Please Note**: $\mu_2$ used in the paper is updated in the FELCEMaterial444.cpp code with the variable mu_neo. It is usually taken as some factor of $\mu_2$ as $\mu_1 = \eta \mu_2$.
