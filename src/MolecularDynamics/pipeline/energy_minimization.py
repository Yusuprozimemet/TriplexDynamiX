from sys import stdout
from openmm import Simulation, StateDataReporter, DCDReporter, PDBReporter

class EnergyMinimization:
    def __init__(self, simulation, config):
        self.simulation = simulation
        self.config = config

    def run_simulation(self):
        # Minimize energy
        self.simulation.reporters.append(StateDataReporter(stdout, self.config['report_interval'], 
                                                           step=True, potentialEnergy=True, temperature=True))
        self.simulation.minimizeEnergy(maxIterations=self.config['max_minimization_iterations'])
        state1 = self.simulation.context.getState(getEnergy=True)
        print(state1.getPotentialEnergy())