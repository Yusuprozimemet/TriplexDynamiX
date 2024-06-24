import yaml

class ConfigLoader:
    def __init__(self, config_file):
        self.config_file = config_file
        self.config = self.load_config()

    def load_config(self):
        with open(self.config_file, 'r') as f:
            config = yaml.safe_load(f)
        return config

    def get_config(self):
        return self.config
