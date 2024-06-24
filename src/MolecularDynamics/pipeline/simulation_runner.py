from sys import stdout
from openmm import Simulation, StateDataReporter, DCDReporter, PDBReporter

class SimulationRunner:
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

        # Run production simulation
        self.simulation.reporters = []
        self.simulation.reporters.append(DCDReporter(self.config['output_dcd_file'], self.config['report_interval']))
        self.simulation.reporters.append(PDBReporter(self.config['output_pdb_file'], self.config['report_interval']))
        self.simulation.reporters.append(StateDataReporter(stdout, self.config['report_interval'], 
                                                           step=True, temperature=True, elapsedTime=True))
        self.simulation.reporters.append(StateDataReporter(self.config['output_csv_file'], self.config['report_interval'], 
                                                           step=True, time=True, potentialEnergy=True, 
                                                           totalEnergy=True, temperature=True))
        self.simulation.step(self.config['simulation_steps'])
