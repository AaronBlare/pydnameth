from pydnameth.model.strategy.load import BetasLoadStrategy
from pydnameth.model.strategy.load import BetasAdjLoadStrategy
from pydnameth.model.strategy.load import BetasHorvathCalculatorLoadStrategy
from pydnameth.model.strategy.load import ObservablesLoadStrategy
from pydnameth.model.strategy.load import ResidualsCommonLoadStrategy
from pydnameth.model.strategy.load import ResidualsSpecialLoadStrategy
from pydnameth.model.strategy.load import EpimutationsLoadStrategy
from pydnameth.model.strategy.load import EntropyLoadStrategy
from pydnameth.model.strategy.load import GenesLoadStrategy
from pydnameth.model.strategy.get import BetasGetStrategy
from pydnameth.model.strategy.get import BetasAdjGetStrategy
from pydnameth.model.strategy.get import BetasHorvathCalculatorGetStrategy
from pydnameth.model.strategy.get import ObservablesGetStrategy
from pydnameth.model.strategy.get import ResidualsCommonGetStrategy
from pydnameth.model.strategy.get import ResidualsSpecialGetStrategy
from pydnameth.model.strategy.get import EpimutationsGetStrategy
from pydnameth.model.strategy.get import EntropyGetStrategy
from pydnameth.model.strategy.get import GenesGetStrategy
from pydnameth.model.strategy.setup import TableSetUpStrategy
from pydnameth.model.strategy.setup import ClockSetUpStrategy
from pydnameth.model.strategy.setup import PlotSetUpStrategy
from pydnameth.model.strategy.setup import CreateSetUpStrategy
from pydnameth.model.strategy.run import TableRunStrategy
from pydnameth.model.strategy.run import ClockRunStrategy
from pydnameth.model.strategy.run import PlotRunStrategy
from pydnameth.model.strategy.run import CreateRunStrategy
from pydnameth.model.strategy.release import ClockReleaseStrategy
from pydnameth.model.strategy.release import TableReleaseStrategy
from pydnameth.model.strategy.release import PlotReleaseStrategy
from pydnameth.model.strategy.release import CreateReleaseStrategy
from pydnameth.model.strategy.save import TableSaveStrategy
from pydnameth.model.strategy.save import ClockSaveStrategy
from pydnameth.model.strategy.save import PlotSaveStrategy
from pydnameth.model.strategy.save import CreateSaveStrategy
from pydnameth.config.experiment.types import Task
from pydnameth.config.experiment.types import DataType


class Context:

    def __init__(self, config):

        if config.experiment.data == DataType.betas:
            self.load_strategy = BetasLoadStrategy()
        elif config.experiment.data == DataType.betas_adj:
            self.load_strategy = BetasAdjLoadStrategy()
        elif config.experiment.data == DataType.betas_horvath_calculator:
            self.load_strategy = BetasHorvathCalculatorLoadStrategy()
        elif config.experiment.data == DataType.observables:
            self.load_strategy = ObservablesLoadStrategy()
        elif config.experiment.data == DataType.residuals_common:
            self.load_strategy = ResidualsCommonLoadStrategy()
        elif config.experiment.data == DataType.residuals_special:
            self.load_strategy = ResidualsSpecialLoadStrategy()
        elif config.experiment.data == DataType.epimutations:
            self.load_strategy = EpimutationsLoadStrategy()
        elif config.experiment.data == DataType.entropy:
            self.load_strategy = EntropyLoadStrategy()
        elif config.experiment.data == DataType.genes:
            self.load_strategy = GenesLoadStrategy()

        if config.experiment.data == DataType.betas:
            self.get_strategy = BetasGetStrategy()
        elif config.experiment.data == DataType.betas_adj:
            self.get_strategy = BetasAdjGetStrategy()
        elif config.experiment.data == DataType.betas_horvath_calculator:
            self.get_strategy = BetasHorvathCalculatorGetStrategy()
        elif config.experiment.data == DataType.observables:
            self.get_strategy = ObservablesGetStrategy()
        elif config.experiment.data == DataType.residuals_common:
            self.get_strategy = ResidualsCommonGetStrategy()
        elif config.experiment.data == DataType.residuals_special:
            self.get_strategy = ResidualsSpecialGetStrategy()
        elif config.experiment.data == DataType.epimutations:
            self.get_strategy = EpimutationsGetStrategy()
        elif config.experiment.data == DataType.entropy:
            self.get_strategy = EntropyGetStrategy()
        elif config.experiment.data == DataType.genes:
            self.get_strategy = GenesGetStrategy()

        if config.experiment.task == Task.table:
            self.setup_strategy = TableSetUpStrategy(self.get_strategy)
        elif config.experiment.task == Task.clock:
            self.setup_strategy = ClockSetUpStrategy(self.get_strategy)
        elif config.experiment.task == Task.plot:
            self.setup_strategy = PlotSetUpStrategy(self.get_strategy)
        elif config.experiment.task == Task.create:
            self.setup_strategy = CreateSetUpStrategy(self.get_strategy)

        if config.experiment.task == Task.table:
            self.run_strategy = TableRunStrategy(self.get_strategy)
        elif config.experiment.task == Task.clock:
            self.run_strategy = ClockRunStrategy(self.get_strategy)
        elif config.experiment.task == Task.plot:
            self.run_strategy = PlotRunStrategy(self.get_strategy)
        elif config.experiment.task == Task.create:
            self.run_strategy = CreateRunStrategy(self.get_strategy)

        if config.experiment.task == Task.table:
            self.release_strategy = TableReleaseStrategy()
        elif config.experiment.task == Task.clock:
            self.release_strategy = ClockReleaseStrategy()
        elif config.experiment.task == Task.plot:
            self.release_strategy = PlotReleaseStrategy()
        elif config.experiment.task == Task.create:
            self.release_strategy = CreateReleaseStrategy()

        if config.experiment.task == Task.table:
            self.save_strategy = TableSaveStrategy()
        elif config.experiment.task == Task.clock:
            self.save_strategy = ClockSaveStrategy()
        elif config.experiment.task == Task.plot:
            self.save_strategy = PlotSaveStrategy()
        elif config.experiment.task == Task.create:
            self.save_strategy = CreateSaveStrategy()

    def pipeline(self, config, configs_child):

        if config.is_run:

            if not self.save_strategy.is_result_exist(config, configs_child):

                config.initialize()
                for config_child in configs_child:
                    config_child.initialize()

                self.load_strategy.load(config, configs_child)
                self.setup_strategy.setup(config, configs_child)
                self.run_strategy.run(config, configs_child)
                self.release_strategy.release(config, configs_child)
                self.save_strategy.save(config, configs_child)
