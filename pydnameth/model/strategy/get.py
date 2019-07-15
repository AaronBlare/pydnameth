import abc
import numpy as np
import pandas as pd
from statsmodels import api as sm


class GetStrategy(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def get_single_base(self, config, items):
        pass

    @abc.abstractmethod
    def get_aux(self, config, item):
        pass

    def get_target(self, config):
        target = config.attributes_dict[config.attributes.target]
        return target


class BetasGetStrategy(GetStrategy):

    def get_single_base(self, config, items):
        rows = [config.base_dict[item] for item in items]
        return config.base_data[np.ix_(rows, config.attributes_indexes)]

    def get_aux(self, config, item):
        aux = ''
        if item in config.cpg_gene_dict:
            aux = ';'.join(config.cpg_gene_dict[item])
        return aux


class BetasAdjGetStrategy(GetStrategy):

    def get_single_base(self, config, items):
        return BetasGetStrategy.get_single_base(self, config, items)

    def get_aux(self, config, item):
        return BetasGetStrategy.get_aux(self, config, item)


class BetasHorvathCalculatorGetStrategy(GetStrategy):

    def get_single_base(self, config, items):
        pass

    def get_aux(self, config, item):
        pass


class BetasSpecGetStrategy(GetStrategy):

    def get_single_base(self, config, items):
        pass

    def get_aux(self, config, item):
        pass


class ResidualsCommonGetStrategy(GetStrategy):

    def get_single_base(self, config, items):
        rows = [config.base_dict[item] for item in items]
        return config.base_data[np.ix_(rows, config.attributes_indexes)]

    def get_aux(self, config, item):
        return BetasGetStrategy.get_aux(self, config, item)


class ResidualsSpecialGetStrategy(GetStrategy):

    def get_single_base(self, config, items):
        cells_dict = config.cells_dict
        exog_df = pd.DataFrame(cells_dict)

        result = np.zeros((len(items), len(config.attributes_indexes)), dtype=np.float32)
        for item_id in range(0, len(items)):
            item = items[item_id]
            row = config.base_dict[item]
            betas = config.base_data[row, config.attributes_indexes]
            endog_dict = {item: betas}
            endog_df = pd.DataFrame(endog_dict)

            reg_res = sm.OLS(endog=endog_df, exog=exog_df).fit()

            residuals = list(map(np.float32, reg_res.resid))

            result[item_id] = residuals

        return result

    def get_aux(self, config, item):
        return BetasGetStrategy.get_aux(self, config, item)


class EpimutationsGetStrategy(GetStrategy):

    def get_single_base(self, config, items):
        rows = [config.betas_dict[item] for item in config.cpg_list if item in config.betas_dict]
        indexes = config.attributes_indexes
        data = np.zeros(len(indexes), dtype=int)

        for subj_id in range(0, len(indexes)):
            col_id = indexes[subj_id]
            subj_col = config.base_data[np.ix_(rows, [col_id])]
            data[subj_id] = np.sum(subj_col)

        data = np.log10(data)
        return [data]

    def get_aux(self, config, item):
        pass


class EntropyGetStrategy(GetStrategy):

    def get_single_base(self, config, items):
        indexes = config.attributes_indexes
        data = config.base_data[indexes]
        return [data]

    def get_aux(self, config, item):
        pass


class ObservablesGetStrategy(GetStrategy):

    def get_single_base(self, config, items):
        pass

    def get_aux(self, config, item):
        pass


class CellsGetStrategy(GetStrategy):

    def get_single_base(self, config, items):
        data = []
        for item in items:
            data.append(config.cells_dict[item])
        return data

    def get_aux(self, config, item):
        pass
