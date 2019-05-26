# -*- coding: utf-8 -*-

"""Top-level package for pydnameth."""
# flake8: noqa

__author__ = """Aaron Blare"""
__email__ = 'aaron.blare@mail.ru'
__version__ = '0.2.3'

from .config.config import Config
from .config.common import CommonTypes
from .config.annotations.annotations import Annotations
from .config.annotations.types import AnnotationKey
from .config.annotations.types import Exclude, CrossReactive, SNP, Chromosome, GeneRegion, Geo, ProbeClass
from .config.attributes.attributes import Cells, Observables, Attributes
from .config.data.data import Data
from .config.data.types import DataPath, DataBase
from .config.experiment.experiment import Experiment
from .config.experiment.types import DataType, Task, Method

from pydnameth.scripts.develop.betas.table import \
    betas_table_aggregator_linreg,\
    betas_table_aggregator_variance,\
    betas_table_linreg,\
    betas_table_cluster
from pydnameth.scripts.develop.betas.clock import \
    betas_clock_linreg,\
    betas_clock_special
from pydnameth.scripts.develop.betas.plot import \
    betas_plot_scatter,\
    betas_plot_curve_clock,\
    betas_plot_variance_histogram_dev

from pydnameth.scripts.develop.betas_horvath_calculator.create import \
    betas_horvath_calculator_create_regular

from pydnameth.scripts.develop.epimutations.load import \
    epimutations_load
from pydnameth.scripts.develop.epimutations.plot import \
    epimutations_plot_scatter,\
    epimutations_plot_range

from pydnameth.scripts.develop.entropy.plot import \
    entropy_plot_scatter

from pydnameth.scripts.develop.genes.table import \
    genes_table_aggregator_linreg,\
    genes_table_aggregator_variance,\
    genes_table_linreg,\
    genes_table_cluster
from pydnameth.scripts.develop.genes.clock import \
    genes_clock_linreg,\
    genes_clock_special
from pydnameth.scripts.develop.genes.plot import \
    genes_plot_scatter,\
    genes_plot_curve_clock,\
    genes_plot_variance_histogram_dev

from pydnameth.scripts.develop.observables.plot import \
    observables_plot_histogram

from pydnameth.scripts.develop.residuals_common.plot import \
    residuals_common_plot_scatter
from pydnameth.scripts.develop.residuals_common.table import \
    residuals_common_table_aggregator_linreg,\
    residuals_common_table_aggregator_variance

from pydnameth.scripts.develop.residuals_special.plot import \
    residuals_special_plot_scatter
from pydnameth.scripts.develop.residuals_special.table import \
    residuals_special_table_aggregator

