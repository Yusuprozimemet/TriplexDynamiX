import os
from datetime import datetime
from openmm import unit, LangevinIntegrator
from openmm.app import HBonds, NoCutoff, PDBFile, Modeller, ForceField, Simulation, StateDataReporter, DCDReporter, PDBReporter
import sys
from molecular_dynamics.utils import find_pdb_file, construct_output_paths

class MolecularDynamicsSimulation:
    def __init__(self, config):
        self.input_folder = config['input_folder']
        self.output_folder = config['output_folder']
        self.forcefield_file = config['forcefield']
        self.temperature = config['temperature']
        self.friction = config['friction']
        self.timestep = config['timestep']
        self.simulation_steps = config['simulation_steps']

        self.input_pdb_path = find_pdb_file(self.input_folder)
        self.base_name = os.path.splitext(os.path.basename(self.input_pdb_path))[0]
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_paths = construct_output_paths(self.output_folder, self.base_name, self.timestamp)

    def run(self):
        pdb = PDBFile(self.input_pdb_path)
        print(pdb.topology)

        modeller = Modeller(pdb.topology, pdb.positions)
        forcefield = ForceField(self.forcefield_file)
        modeller.addHydrogens(forcefield)

        with open(self.output_paths['output_init_with_hydrogens_path'], 'w') as outfile:
            PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

        print(f"Modified PDB file with hydrogens saved to: {self.output_paths['output_init_with_hydrogens_path']}")

        system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=HBonds)

        integrator = LangevinIntegrator(self.temperature*unit.kelvin, self.friction/unit.picosecond, self.timestep*unit.femtoseconds)

        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

        simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))
        simulation.minimizeEnergy(maxIterations=1000)
        state1 = simulation.context.getState(getEnergy=True)
        print(state1.getPotentialEnergy())

        simulation.reporters = []

        simulation.reporters.append(DCDReporter(self.output_paths['output_dcd_path'], 1000))
        simulation.reporters.append(PDBReporter(self.output_paths['output_pdb_path'], 1000))
        simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, temperature=True, elapsedTime=True))
        simulation.reporters.append(StateDataReporter(self.output_paths['output_csv_path'], 1000, step=True, time=True,
                                                      potentialEnergy=True, totalEnergy=True, temperature=True))

        simulation.step(self.simulation_steps)

        del simulation

        print(f"Output PDB file saved to: {self.output_paths['output_pdb_path']}")
        print(f"Output DCD file saved to: {self.output_paths['output_dcd_path']}")
        print(f"Output CSV file saved to: {self.output_paths['output_csv_path']}")
