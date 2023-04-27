import awkward as ak
import copy

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.parameters.lumi import goldenJSON
from pocket_coffea.parameters.btag import btag

class genWeightsProcessor(BaseProcessorABC):
    def __init__(self, cfg: Configurator) -> None:
        super().__init__(cfg)
        self.output_format = {
            "sum_genweights": {},
            "cutflow": {
                "initial": {s: 0 for s in self.cfg.datasets}
            }
        }

    def load_metadata(self):
        self._dataset = self.events.metadata["dataset"]
        self._sample = self.events.metadata["sample"]
        self._year = self.events.metadata["year"]
        self._btag = btag[self._year]
        self._isMC = self.events.metadata["isMC"] == "True"
        
        if self._isMC:
            self._era = "MC"
            self._xsec = self.events.metadata["xsec"]
        else:
            self._era = self.events.metadata["era"]
            self._goldenJSON = goldenJSON[self._year]

    def apply_object_preselection(self, variation):
        pass

    def count_objects(self, variation):
        pass

    def process(self, events: ak.Array):
        self.events = events

        self.load_metadata()
        self.output = copy.deepcopy(self.output_format)

        self.nEvents_initial = self.nevents
        self.output['cutflow']['initial'][self._dataset] += self.nEvents_initial
        if self._isMC:
            self.output['sum_genweights'][self._dataset] = ak.sum(self.events.genWeight)

        return self.output
