from project_name.config_loader import ConfigLoader
from project_name.simulation_loader import PDBForceFieldLoader, SimulationSetup
from project_name.simulation_runner import SimulationRunner

def main():
    # Load config
    config_loader = ConfigLoader('config/config.yaml')
    config = config_loader.get_config()

    # Load PDB and force field
    pdb_forcefield_loader = PDBForceFieldLoader(config)
    modeller, forcefield = pdb_forcefield_loader.modeller, pdb_forcefield_loader.forcefield

    # Setup simulation
    simulation_setup = SimulationSetup(config, modeller, forcefield)
    simulation = simulation_setup.setup_simulation()

    # Run simulation
    simulation_runner = SimulationRunner(simulation, config)
    simulation_runner.run_simulation()

if __name__ == "__main__":
    main()
