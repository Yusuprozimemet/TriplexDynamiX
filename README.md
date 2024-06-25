# MolecularDynamicsApp

Welcome to our molecular dynamics application designed for trajectory analysis, energy minimization, and comprehensive data handling and interpretation. More features are currently in development. Initially, this app was used for analyzing RNA triples stability. Although not all parameters and features are included at present, it serves as a useful starting point for researchers new to molecular dynamics. For further applications, refer to the OpenMM and MDtraj documentation.

### How to Use This App:

1. Clone the repository:
   ```
   git clone https://github.com/Yusuprozimemet/MolecularDynamicsApp.git
   ```

2. Create a Python 3.7 environment:
   ```
   conda create -n yourenv python=3.7
   ```

3. Activate your environment:
   ```
   conda activate yourenv
   ```

4. Update your environment:
   ```
   conda env update -n yourenv -f environment.yaml
   ```

5. Place your input PDB files in `src/MolecularDynamics/data/input`. Ensure these files are ready for molecular dynamics simulations. You will find the output in `src/MolecularDynamics/data/output`.

6. Navigate to the appropriate directory to run the desired Python script.

7. You can modify parameters in the YAML files for the scripts. 

8. For any questions or issues, please contact me at Rouzimaimaiti.Yusufu@mail.huji.ac.il.

Enjoy! 