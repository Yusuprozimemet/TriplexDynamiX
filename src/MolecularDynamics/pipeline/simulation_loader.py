import os
from openmm import app, unit
from openmm.app import PDBFile, Modeller
from openmm.app.forcefield import ForceField
from openmm.unit import *

class PDBForceFieldLoader:
    def __init__(self, config):
        self.config = config
        self.modeller, self.forcefield = self.load_pdb_and_forcefield()

    def load_pdb_and_forcefield(self):
        pdb_path = self.config['pdb_path']
        if not os.path.isfile(pdb_path):
            raise FileNotFoundError(f"File not found: {pdb_path}")

        # Load PDB file and print topology
        pdb = PDBFile(pdb_path)
        modeller = Modeller(pdb.topology, pdb.positions)

        # Load force field
        forcefield = ForceField(*self.config['forcefield_files'])
        modeller.addHydrogens(forcefield)

        return modeller, forcefield

class SimulationSetup:
    def __init__(self, config, modeller, forcefield):
        self.config = config
        self.modeller = modeller
        self.forcefield = forcefield

    def setup_simulation(self):
        system = self.forcefield.createSystem(self.modeller.topology, 
                                              nonbondedMethod=app.NoCutoff, 
                                              nonbondedCutoff=3*nanometer, 
                                              constraints=app.HBonds)
        integrator = LangevinIntegrator(self.config['temperature']*kelvin, 
                                        self.config['friction_coefficient']/picosecond, 
                                        self.config['time_step']*femtoseconds)
        simulation = Simulation(self.modeller.topology, system, integrator)
        simulation.context.setPositions(self.modeller.positions)
        return simulation
