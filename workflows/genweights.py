import copy

import awkward as ak

from workflows.fatjet_base import fatjetBaseProcessor

class genWeightsProcessor(fatjetBaseProcessor):

    def process(self, events: ak.Array):
        self.events = events

        ###################
        # At the beginning of the processing the initial number of events
        # and the sum of the genweights is stored for later use
        #################
        self.load_metadata()
        # Define the accumulator instance for this chunk
        self.output = copy.deepcopy(self.output_format)

        self.nEvents_initial = self.nevents
        self.output['cutflow']['initial'][self._sample] += self.nEvents_initial
        if self._isMC:
            # This is computed before any preselection
            key = f"{self._sample}_{self._part}"
            self.output['sum_genweights'][key] = ak.sum(self.events.genWeight)
            if self._hasSubsamples:
                for subs in self._subsamples_names:
                    self.output['sum_genweights'][subs] = self.output['sum_genweights'][key]

        return self.output
