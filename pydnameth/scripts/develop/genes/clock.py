import copy
import os
from shutil import copyfile
from anytree import Node
from pydnameth.config.config import Config
from pydnameth.config.experiment.types import Task, Method, DataType
from pydnameth.config.experiment.experiment import Experiment
from pydnameth.infrastucture.file_name import get_file_name
from pydnameth.infrastucture.path import get_save_path
from pydnameth.model.tree import build_tree, calc_tree


def genes_clock_linreg(
    data,
    annotations,
    attributes,
    method_params=None
):
    """
        Temporary description
    """

    config_root = Config(
        data=copy.deepcopy(data),
        experiment=Experiment(
            data=DataType.genes,
            task=Task.clock,
            method=Method.linreg,
            method_params=copy.deepcopy(method_params)
        ),
        annotations=copy.deepcopy(annotations),
        attributes=copy.deepcopy(attributes),
        is_run=True,
        is_root=True
    )
    root = Node(name=str(config_root), config=config_root)

    config_child = Config(
        data=copy.deepcopy(data),
        experiment=Experiment(
            data=DataType.genes,
            task=Task.table,
            method=Method.linreg
        ),
        annotations=copy.deepcopy(annotations),
        attributes=copy.deepcopy(attributes),
        is_run=True,
        is_root=False
    )

    Node(name=str(config_child), config=config_child, parent=root)

    build_tree(root)
    calc_tree(root)


def genes_clock_special(
    data,
    annotations,
    attributes,
    file,
    method_params=None,
):
    """
        Temporary description
    """

    if os.path.isfile(file):

        head, tail = os.path.split(file)
        fn = os.path.splitext(tail)[0]
        ext = os.path.splitext(tail)[1]

        config_root = Config(
            data=copy.deepcopy(data),
            experiment=Experiment(
                data=DataType.genes,
                task=Task.clock,
                method=Method.linreg,
                method_params=copy.deepcopy(method_params)
            ),
            annotations=copy.deepcopy(annotations),
            attributes=copy.deepcopy(attributes),
            is_run=True,
            is_root=True
        )
        root = Node(name=str(config_root), config=config_root)

        config_child = Config(
            data=copy.deepcopy(data),
            experiment=Experiment(
                data=DataType.genes,
                task=Task.table,
                method=Method.special,
                method_params={'file_name': fn}
            ),
            annotations=copy.deepcopy(annotations),
            attributes=copy.deepcopy(attributes),
            is_run=False,
            is_root=False
        )

        Node(name=str(config_child), config=config_child, parent=root)

        build_tree(root)

        new_file = get_save_path(config_child) + '/' + \
            get_file_name(config_child) + ext

        copyfile(file, new_file)

        calc_tree(root)

    else:
        raise FileNotFoundError(f'File {file} not found.')
