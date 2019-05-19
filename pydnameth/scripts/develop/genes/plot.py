import copy
from anytree import Node
from pydnameth.config.config import Config
from pydnameth.config.experiment.types import Task, Method, DataType
from pydnameth.config.experiment.experiment import Experiment
from pydnameth.config.attributes.attributes import Observables, Cells, Attributes
from pydnameth.model.tree import build_tree, calc_tree
from pydnameth.scripts.develop.plot import plot_scatter


def genes_plot_scatter(
    data,
    annotations,
    attributes,
    genes_list,
    observables_list,
    child_method=Method.linreg,
    method_params=None
):
    """
    Plotting methylation level from observables as scatter for provided subjects subsets and provided CpG list.

    Possible parameters of experiment:

     * ``'x_range'``: can be ``'auto'`` or list with two elements, which are borders of target axis.
     * ...

    :param data: pdm.Data instance, which specifies information about dataset.
    :param annotations: pdm.Annotations instance, which specifies subset of CpGs.
    :param attributes: pdm.Attributes instance, which specifies information about subjects.
    :param cpg_list: List of CpGs for plotting
    :param observables_list: list of subjects subsets. Each element in list is dict,
     where ``key`` is observable name and ``value`` is possible values for this observable.
    :param method_params: parameters of experiment.
    """

    plot_scatter(
        data_type=DataType.genes,
        data=data,
        annotations=annotations,
        attributes=attributes,
        cpg_list=genes_list,
        observables_list=observables_list,
        child_method=child_method,
        method_params=method_params
    )


def genes_plot_curve_clock(
    data,
    annotations,
    attributes,
    observables_list,
    child_method=Method.linreg,
    data_params=None,
    method_params=None
):
    data_type = DataType.genes

    clock_method_params = {
        'type': 'all',
        'part': 0.25,
        'size': 100,
        'runs': 100,
    }

    config_root = Config(
        data=copy.deepcopy(data),
        experiment=Experiment(
            data=data_type,
            task=Task.plot,
            method=Method.curve,
            data_params=copy.deepcopy(data_params),
            method_params=copy.deepcopy(method_params)
        ),
        annotations=copy.deepcopy(annotations),
        attributes=copy.deepcopy(attributes),
        is_run=True,
        is_root=True,
        is_load_child=True
    )

    root = Node(name=str(config_root), config=config_root)

    for d in observables_list:

        observables_child = Observables(
            name=copy.deepcopy(attributes.observables.name),
            types=d
        )

        cells_child = Cells(
            name=copy.deepcopy(attributes.cells.name),
            types=copy.deepcopy(attributes.cells.types)
        )

        attributes_child = Attributes(
            target=copy.deepcopy(attributes.target),
            observables=observables_child,
            cells=cells_child,
        )

        config_child_lvl_1 = Config(
            data=copy.deepcopy(data),
            experiment=Experiment(
                data=data_type,
                task=Task.clock,
                method=copy.deepcopy(child_method),
                data_params=copy.deepcopy(data_params),
                method_params=clock_method_params
            ),
            annotations=copy.deepcopy(annotations),
            attributes=attributes_child,
            is_run=False,
            is_root=False,
            is_load_child=False
        )
        node_lvl_1 = Node(name=str(config_child_lvl_1), config=config_child_lvl_1, parent=root)

        config_child_lvl_2 = Config(
            data=copy.deepcopy(data),
            experiment=Experiment(
                data=data_type,
                task=Task.table,
                method=Method.linreg
            ),
            annotations=copy.deepcopy(annotations),
            attributes=attributes_child,
            is_run=False,
            is_root=False
        )
        Node(name=str(config_child_lvl_2), config=config_child_lvl_2, parent=node_lvl_1)

    build_tree(root)
    calc_tree(root)


def genes_plot_variance_histogram_dev(
    data,
    annotations,
    attributes,
    genes_list,
    observables_list,
    child_method=Method.linreg,
    method_params=None
):
    for gene in genes_list:

        config_root = Config(
            data=copy.deepcopy(data),
            experiment=Experiment(
                data=DataType.genes,
                task=Task.methylation,
                method=Method.variance_histogram,
                method_params=copy.deepcopy(method_params)
            ),
            annotations=copy.deepcopy(annotations),
            attributes=copy.deepcopy(attributes),
            is_run=True,
            is_root=True,
            is_load_child=False
        )

        if config_root.experiment.method_params is None:
            config_root.experiment.method_params = dict()

        config_root.experiment.method_params['item'] = gene

        root = Node(name=str(config_root), config=config_root)

        for d in observables_list:
            observables_child = Observables(
                name=copy.deepcopy(attributes.observables.name),
                types=d
            )

            cells_child = Cells(
                name=copy.deepcopy(attributes.cells.name),
                types=copy.deepcopy(attributes.cells.types)
            )

            attributes_child = Attributes(
                target=copy.deepcopy(attributes.target),
                observables=observables_child,
                cells=cells_child,
            )

            config_child = Config(
                data=copy.deepcopy(data),
                experiment=Experiment(
                    data=DataType.genes,
                    task=Task.table,
                    method=copy.deepcopy(child_method)
                ),
                annotations=copy.deepcopy(annotations),
                attributes=attributes_child,
                is_run=False,
                is_root=False,
                is_load_child=False
            )
            Node(name=str(config_child), config=config_child, parent=root)

        build_tree(root)
        calc_tree(root)
