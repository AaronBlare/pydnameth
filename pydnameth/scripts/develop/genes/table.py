from pydnameth.config.experiment.types import Method, DataType
from pydnameth.scripts.develop.table import table, table_aggregator_linreg, table_aggregator_variance


def genes_table_linreg(
    data,
    annotations,
    attributes,
    method_params=None
):
    """
        Temporary description
    """

    table(
        data=data,
        annotations=annotations,
        attributes=attributes,
        data_type=DataType.genes,
        method=Method.linreg,
        method_params=method_params,
    )


def genes_table_cluster(
    data,
    annotations,
    attributes,
    method_params=None
):
    table(
        data=data,
        annotations=annotations,
        attributes=attributes,
        data_type=DataType.genes,
        method=Method.cluster,
        method_params=method_params,
    )


def genes_table_aggregator_linreg(
    data,
    annotations,
    attributes,
    observables_list,
    method_params=None
):
    """
        Temporary description
    """

    table_aggregator_linreg(
        DataType.genes,
        data,
        annotations,
        attributes,
        observables_list,
        method_params
    )


def genes_table_aggregator_variance(
    data,
    annotations,
    attributes,
    observables_list,
    data_params,
):
    """
        Temporary description
    """

    table_aggregator_variance(
        DataType.genes,
        data,
        annotations,
        attributes,
        observables_list,
        data_params=data_params,
    )
