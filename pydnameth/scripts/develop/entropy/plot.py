import copy
from anytree import Node
from pydnameth.config.config import Config
from pydnameth.config.experiment.types import Task, Method, DataType
from pydnameth.config.experiment.experiment import Experiment
from pydnameth.config.attributes.attributes import Observables, Cells, Attributes
from pydnameth.model.tree import build_tree, calc_tree


def entropy_plot_scatter(
    data,
    annotations,
    attributes,
    observables_list,
    data_params=None,
    method_params=None
):
    entropy_plot(
        data=data,
        annotations=annotations,
        attributes=attributes,
        observables_list=observables_list,
        method=Method.scatter,
        data_params=data_params,
        method_params=method_params
    )


def entropy_plot(
    data,
    annotations,
    attributes,
    observables_list,
    method,
    data_params=None,
    method_params=None
):
    config_root = Config(
        data=copy.deepcopy(data),
        experiment=Experiment(
            data=DataType.entropy,
            task=Task.plot,
            method=method,
            data_params=copy.deepcopy(data_params),
            method_params=copy.deepcopy(method_params)
        ),
        annotations=copy.deepcopy(annotations),
        attributes=copy.deepcopy(attributes),
        is_run=True,
        is_root=True,
        is_load_child=False
    )

    root = Node(name=str(config_root), config=config_root)

    for types in observables_list:
        observables_child = Observables(
            name=copy.deepcopy(attributes.observables.name),
            types=types
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
                data=DataType.entropy,
                task=Task.table,
                method=Method.mock,
                data_params=copy.deepcopy(data_params)
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
