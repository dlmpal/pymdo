from typing import Dict
from typing import List
from copy import deepcopy
from os.path import exists

from numpy import ndarray, array
import json
import h5py

from pymdo.core.variable import Variable

from .variable import Variable
from .variable import check_values_match
from .variable import FLOAT_DATA_TYPE


class Cache:

    CACHE_POLICY_LATEST = "latest"

    CACHE_POLICY_FULL = "full"

    CACHE_POLICY = [CACHE_POLICY_LATEST,
                    CACHE_POLICY_FULL]

    def __init__(self,
                 inputVars: List[Variable],
                 outputVars: List[Variable],
                 policy: str = CACHE_POLICY_LATEST,
                 tol: float = 1e-9,
                 path: str = None) -> None:

        self.inputVars = inputVars

        self.outputVars = outputVars

        self.diffInputs = [var for var in self.inputVars]

        self.diffOutputs = [var for var in self.outputVars]

        self.policy = policy

        self.tol = tol

        self.path = path

        self.from_disk()
    
    def update_diff_inputs(self, vars: List[Variable]) -> None:

        self.diffInputs = [var for var in vars]
    
    def update_diff_outputs(self, vars: List[Variable]) -> None:

        self.diffOutputs = [var for var in vars]

    def check_if_entry_exists(self, inputs: Dict[str, ndarray]):
        """
        Check if there exists an entry for the given inputs

        """

        raise NotImplementedError

    def add_entry(self,
                  inputs: Dict[str, ndarray],
                  outputs: Dict[str, ndarray] = None,
                  jac: Dict[str, Dict[str, ndarray]] = None) -> None:
        """

        Adds an entry in the cache for the given inputs.

        If an entry already exists, the existing output/jacobian values 
        are overriden by the ones provided.

        """

        raise NotImplementedError

    def load_entry(self, inputs: Dict[str, ndarray]):
        """

        Load the entry for the given inputs, if it exists

        """

        raise NotImplementedError

    def from_disk(self, path: str = None):
        """

        Load the cache entries from disk

        """

        raise NotImplementedError

    def to_disk(self, path: str = None):
        """

        Save the cache entries to disk

        """

        raise NotImplementedError

class MemoryCache(Cache):

    def __init__(self,
                 inputVars: List[Variable],
                 outputVars: List[Variable],
                 policy: str = Cache.CACHE_POLICY_LATEST,
                 tol: float = 1e-9,
                 path: str = None) -> None:
        
        self.entries = []

        super().__init__(inputVars,
                         outputVars,
                         policy,
                         tol,
                         path)

    def check_if_entry_exists(self, inputs: Dict[str, ndarray]):

        if not self.entries:

            return None

        for entry in reversed(self.entries):

            if check_values_match(self.inputVars,
                                  entry["inputs"],
                                  inputs,
                                  self.tol):

                return entry

        return None

    def add_entry(self,
                  inputs: Dict[str, ndarray],
                  outputs: Dict[str, ndarray] = None,
                  jac: Dict[str, Dict[str, ndarray]] = None) -> None:
        """

        Adds an entry in the cache for the given inputs.

        If an entry already exists, the existing output/jacobian values 
        are overriden by the ones provided.

        """

        newEntry = self.check_if_entry_exists(inputs)

        entryExists = True

        # If entry does not exist, instantiate new one
        if newEntry is None:

            newEntry = {"inputs": {var.name: deepcopy(inputs[var.name]) for var in self.inputVars},
                        "outputs": {},
                        "jac": {}}

            entryExists = False

        # Modify entry outputs, if given
        if outputs is not None:
            newEntry["outputs"] = {var.name: deepcopy(
                outputs[var.name]) for var in self.outputVars}

        # Modify entry jac, if given
        if jac is not None:

            for outVar in self.diffOutputs:

                if outVar.name in jac:

                    newEntry["jac"][outVar.name] = {}

                for inVar in self.diffInputs:

                    if inVar.name in jac[outVar.name]:

                        newEntry["jac"][outVar.name][inVar.name] = deepcopy(
                            jac[outVar.name][inVar.name])

        if self.policy == self.CACHE_POLICY_LATEST:

            self.entries = []

            entryExists = False

        if entryExists == False:

            self.entries.append(newEntry)

    def load_entry(self, inputs: Dict[str, ndarray]):

        entry = self.check_if_entry_exists(inputs)

        if entry is not None:

            return deepcopy(entry["outputs"]), deepcopy(entry["jac"])

        return None, None

    def from_disk(self, path: str = None):

        if path is None:

            path = self.path
        
        if not exists(path + ".json"):

            return

        with open(path + ".json", "r") as file:

            jsonObj = json.load(file)

            self.entries = []

            for _, entry in jsonObj.items():

                for inVar in entry["inputs"]:

                    entry["inputs"][inVar] = array(entry["inputs"][inVar],
                                                   dtype=FLOAT_DATA_TYPE)

                for outVar in entry["outputs"]:

                    entry["outputs"][outVar] = array(entry["outputs"][outVar],
                                                     dtype=FLOAT_DATA_TYPE)

                for outVarName in entry["jac"]:

                    for invarName in entry["jac"][outVarName]:

                        entry["jac"][outVarName][invarName] = array(entry["jac"][outVarName][invarName],
                                                                    dtype=FLOAT_DATA_TYPE)

                self.entries.append(entry)

    def to_disk(self, path: str = None):

        if path is None:

            path = self.path

        jsonObj = {}

        for idx, entry in enumerate(self.entries):

            jsonObj[idx] = {"inputs": {},
                            "outputs": {},
                            "jac": {}}

            for inVar in self.inputVars:

                jsonObj[idx]["inputs"][inVar.name] = entry["inputs"][inVar.name].tolist()

            for outVar in self.outputVars:

                jsonObj[idx]["outputs"][outVar.name] = entry["outputs"][outVar.name].tolist()

            for outVarName in entry["jac"]:

                jsonObj[idx]["jac"][outVarName] = {}

                for inVarName in entry["jac"][outVarName]:

                    jsonObj[idx]["jac"][outVarName][inVarName] = entry["jac"][outVarName][inVarName].tolist(
                    )

        jsonObj = json.dumps(jsonObj, indent=4)

        with open(path+".json", "w") as file:

            file.write(jsonObj)

class FileCache(Cache):

    def __init__(self,
                 inputVars: List[Variable],
                 outputVars: List[Variable],
                 policy: str = Cache.CACHE_POLICY_LATEST,
                 tol: float = 1e-9,
                 path: str = None) -> None:
        
        self.entries: h5py.File = None

        if not exists(path + ".hdf5"):

            self.entries = h5py.File(path + ".hdf5", "w")
        
        super().__init__(inputVars,
                         outputVars,
                         policy,
                         tol,
                         path)
        
    def check_if_entry_exists(self, inputs: Dict[str, ndarray]):

        if not self.entries:

            return None

        for entryIdx in reversed(self.entries):

            if check_values_match(self.inputVars,
                                  self.entries[entryIdx]["inputs"],
                                  inputs,
                                  self.tol):

                return self.entries[entryIdx]

        return None

    def add_entry(self,
                  inputs: Dict[str, ndarray],
                  outputs: Dict[str, ndarray] = None,
                  jac: Dict[str, Dict[str, ndarray]] = None) -> None:

        newEntry = None

        if self.policy == Cache.CACHE_POLICY_LATEST:
            if self.entries:
                del self.entries["0"]
        else:
            newEntry = self.check_if_entry_exists(inputs)

        # If entry does not exist, instantiate new one
        if newEntry is None:

            newEntry = self.entries.create_group(str(len(self.entries)))

            newEntry.create_group("inputs")

            newEntry.create_group("outputs")

            newEntry.create_group("jac")

            for var in self.inputVars:

                newEntry["inputs"].create_dataset(
                    var.name, 
                    data=inputs[var.name],
                    dtype = FLOAT_DATA_TYPE)

        # Modify entry outputs, if given
        if outputs is not None:

            for var in self.outputVars:

                newEntry["outputs"].require_dataset(name=var.name,
                                                    data=outputs[var.name],
                                                    shape=(var.size),
                                                    dtype=FLOAT_DATA_TYPE)

        # Modify entry jac, if given
        if jac is not None:

            for outVar in self.diffOutputs:

                if outVar.name in jac:

                    newEntry["jac"].require_group(outVar.name)

                for inVar in self.diffInputs:

                    if inVar.name in jac[outVar.name]:

                        newEntry["jac"][outVar.name].require_dataset(name=inVar.name,
                                                                     data=jac[outVar.name][inVar.name],
                                                                     shape=(
                                                                         outVar.size, inVar.size),
                                                                     dtype=FLOAT_DATA_TYPE)

    def load_entry(self, inputs: Dict[str, ndarray]):

        entry = self.check_if_entry_exists(inputs)

        if entry is not None:

            outputs = None 

            if len(entry["outputs"]) > 0:

                outputs = {var.name: entry["outputs"][var.name][()] for var in self.outputVars}

            jac = None 

            if len(entry["jac"]) > 0:

               jac = {outVarName: {inVarName: entry["jac"][outVarName][inVarName][()] for inVarName in entry["jac"][outVarName].keys()}
                      for outVarName in entry["jac"].keys()}

            return outputs, jac

        return None, None

    def from_disk(self, path: str = None):

        if path is None:

            path = self.path

        if not exists(path + ".hdf5"):
            return
        
        self.entries = h5py.File(path + ".hdf5", "r+")

    def to_disk(self, path: str = None):

        pass

def cache_factory(inputVars: List[Variable],
                  outputVars: List[Variable],
                  type: str = "memory",
                  policy: str = "latest",
                  tol: float = 1e-9,
                  path: str = None):

    if type is None:

        return None

    if type == "memory":

        return MemoryCache(inputVars,
                           outputVars,
                           policy,
                           tol,
                           path)

    elif type == "file":

        return FileCache(inputVars,
                         outputVars,
                         policy,
                         tol,
                         path)
