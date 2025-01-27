import logging


class Logger:
    def __init__(self, level=logging.INFO, silent=True):
        self.logger = logging.getLogger("Shampoo")
        self.logger.setLevel(level)
        self.silent = silent

        if not self.silent:
            console = logging.StreamHandler()
            console.setLevel(level)
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            console.setFormatter(formatter)
            self.logger.addHandler(console)

    def log(self, message, level=logging.INFO):
        if not self.silent:
            self.logger.log(level, message)
