# simulation.py

import os
from datetime import datetime
from openmm import unit, LangevinIntegrator
from openmm.app import HBonds, NoCutoff, PDBFile, Modeller, ForceField, Simulation, StateDataReporter, DCDReporter, PDBReporter
import sys
import csv
from energy_minimization.utils import find_pdb_file, construct_output_paths

class MinimizeEnergy:
    def __init__(self, config, input_pdb_path, output_paths):
        self.input_pdb_path = input_pdb_path
        self.output_paths = output_paths

        self.forcefield_file = config['forcefield']
        self.temperature = config['temperature']
        self.friction = config['friction']
        self.timestep = config['timestep']

        self.base_name = os.path.splitext(os.path.basename(self.input_pdb_path))[0]

    def run(self):
        pdb = PDBFile(self.input_pdb_path)
        print(pdb.topology)

        modeller = Modeller(pdb.topology, pdb.positions)
        forcefield = ForceField(self.forcefield_file)
        modeller.addHydrogens(forcefield)

        system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=HBonds)

        integrator = LangevinIntegrator(self.temperature*unit.kelvin, self.friction/unit.picosecond, self.timestep*unit.femtoseconds)

        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

        print("Minimizing energy...")
        simulation.minimizeEnergy(maxIterations=1000)

        state = simulation.context.getState(getEnergy=True)
        minimized_energy = state.getPotentialEnergy()

        print(f"Minimized Potential Energy: {minimized_energy}")

        # Write minimized energy to CSV file
        csv_file = self.output_paths['output_csv_path']
        with open(csv_file, 'a', newline='') as csvfile:
            fieldnames = ['pdb_file', 'minimized_energy']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            # Write header if the file is empty
            if os.stat(csv_file).st_size == 0:
                writer.writeheader()

            writer.writerow({'pdb_file': os.path.basename(self.input_pdb_path), 'minimized_energy': minimized_energy / unit.kilojoules_per_mole})

        print(f"Minimized energy saved to: {csv_file}")

        del simulation
