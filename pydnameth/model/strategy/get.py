import abc
import numpy as np


class GetStrategy(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def get_single_base(self, config, items):
        pass

    @abc.abstractmethod
    def get_aux(self, config, item):
        pass

    def get_target(self, config, item='any', categorical=True):
        tmp = np.array(self.get_observalbe(config, key=config.attributes.target, item=item, categorical=categorical))
        return tmp

    def get_observalbe(self, config, key, item='any', categorical=True):
        if config.base_missed_dict is not None and len(config.base_missed_dict[item]) > 0:
            passed_ids = []
            for id, col in enumerate(config.attributes_indexes):
                if col not in config.base_missed_dict[item]:
                    passed_ids.append(id)
            data = []
            for id in passed_ids:
                if categorical:
                    data.append(config.observables_categorical_dict[key][id])
                else:
                    data.append(config.observables_dict[key][id])
        else:
            if categorical:
                data = config.observables_categorical_dict[key]
            else:
                data = config.observables_dict[key]
        return data

    def get_cell(self, config, key, item='any'):
        if config.base_missed_dict is not None and len(config.base_missed_dict[item]) > 0:
            passed_ids = []
            for id, col in enumerate(config.attributes_indexes):
                if col not in config.base_missed_dict[item]:
                    passed_ids.append(id)
            data = []
            for id in passed_ids:
                data.append(config.cells_dict[key][id])
        else:
            data = config.cells_dict[key]
        return data


class BetasGetStrategy(GetStrategy):

    def get_single_base(self, config, item):
        row = config.base_dict[item]
        if len(config.base_missed_dict[item]) > 0:
            cols = []
            for col in config.attributes_indexes:
                if col not in config.base_missed_dict[item]:
                    cols.append(col)
            data = config.base_data[row, cols]
        else:
            data = config.base_data[row, config.attributes_indexes]
        return data

    def get_aux(self, config, item):
        aux = ''
        if item in config.cpg_gene_dict:
            aux = ';'.join(config.cpg_gene_dict[item])
        return aux


class BOPsGetStrategy(GetStrategy):

    def get_single_base(self, config, item):
        row = config.target_dict[item]
        if len(config.base_missed_dict[item]) > 0:
            cols = []
            for col in config.attributes_indexes:
                if col not in config.base_missed_dict[item]:
                    cols.append(col)
            data = config.base_data[row, cols]
        else:
            data = config.base_data[row, config.attributes_indexes]
        return data

    def get_aux(self, config, item):
        cpgs = config.bops[item]['cpg']
        aux = ';'.join(cpgs)
        return aux


class BetasAdjGetStrategy(GetStrategy):

    def get_single_base(self, config, item):
        return BetasGetStrategy.get_single_base(self, config, item)

    def get_aux(self, config, item):
        return BetasGetStrategy.get_aux(self, config, item)


class BetasHorvathCalculatorGetStrategy(GetStrategy):

    def get_single_base(self, config, item):
        pass

    def get_aux(self, config, item):
        pass


class BetasSpecGetStrategy(GetStrategy):

    def get_single_base(self, config, item):
        pass

    def get_aux(self, config, item):
        pass


class ResidualsGetStrategy(GetStrategy):

    def get_single_base(self, config, item):
        return BetasGetStrategy.get_single_base(self, config, item)

    def get_aux(self, config, item):
        return BetasGetStrategy.get_aux(self, config, item)


class ResidOldGetStrategy(GetStrategy):

    def get_single_base(self, config, item):
        return BetasGetStrategy.get_single_base(self, config, item)

    def get_aux(self, config, item):
        return BetasGetStrategy.get_aux(self, config, item)


class EpimutationsGetStrategy(GetStrategy):

    def get_single_base(self, config, item):
        rows = [config.betas_dict[cpg] for cpg in config.cpg_list if cpg in config.betas_dict]
        indexes = config.attributes_indexes
        data = np.zeros(len(indexes), dtype=int)

        for subj_id in range(0, len(indexes)):
            col_id = indexes[subj_id]
            subj_col = config.base_data[np.ix_(rows, [col_id])]
            data[subj_id] = np.sum(subj_col)

        data = np.log10(data)
        return data

    def get_aux(self, config, item):
        pass


class EntropyGetStrategy(GetStrategy):

    def get_single_base(self, config, item):
        indexes = config.attributes_indexes
        data = config.base_data[indexes]
        return data

    def get_aux(self, config, item):
        pass


class ObservablesGetStrategy(GetStrategy):

    def get_single_base(self, config, item):
        pass

    def get_aux(self, config, item):
        pass


class CellsGetStrategy(GetStrategy):

    def get_single_base(self, config, item):
        if len(config.base_missed_dict[item]) > 0:
            data = []
            for col in config.attributes_indexes:
                if col not in config.base_missed_dict[item]:
                    data.append(config.cells_dict[item][col])
        else:
            data = config.cells_dict[item]
        return data

    def get_aux(self, config, item):
        pass


class GenesGetStrategy(GetStrategy):

    def get_single_base(self, config, item):
        return BetasGetStrategy.get_single_base(self, config, item)

    def get_aux(self, config, item):
        aux = ''
        return aux
