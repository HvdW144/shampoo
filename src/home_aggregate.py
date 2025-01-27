import warnings

warnings.filterwarnings("ignore")

import numpy as np


# Shampoo HomeAggregate file.
#
# Written and documented by Mark Oosterloo
#
# Version: 15-03-2024

class HomeAggregate:

    def __init__(self, model, size):
        """
        Initializes a variety of parameters associated with the home aggregate.
        """

        self.initProps(size, model)

    def initProps(self, size, model):
        self.prop = {}

        if model.pisoBenchmark:
            self.prop["phi"] = 1
            self.prop["rhoAgg"] = 2000
        elif model.phi != None:
            self.prop["phi"] = model.phi
            self.prop["rhoAgg"] = 2094 * self.prop["phi"]
        else:
            self.prop["phi"] = float(model.paraDict["phi"])
            self.prop["rhoAgg"] = 2094 * self.prop["phi"]
        self.prop["sAgg"] = size  # in m; initial size is equal to the monomer size.

        self.prop["mAgg"] = 4 / 3 * np.pi * self.prop["sAgg"] ** 3 * self.prop["rhoAgg"]
