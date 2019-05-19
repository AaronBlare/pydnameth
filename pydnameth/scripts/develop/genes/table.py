from pydnameth.config.experiment.types import Method, DataType
from pydnameth.scripts.develop.table import table, table_aggregator_linreg, table_aggregator_variance


def genes_table_linreg(
    data,
    annotations,
    attributes,
    method_params=None
):
    """
        There should be a description
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
        Producing table with information about observable-specificity of target data type
        and target observable for each CpG.

        Columns:

        * item: CpG id.
        * aux: gene, on which CpG is mapped.
        * area_intersection_rel: relative intersection area of polygons
          which is equals area of polygon(s) intersection to area of polygons union ratio.
        * slope_intersection_rel: relative intersection area of allowed regions for slopes of linear regression.
        * max_abs_slope: maximal absolute slope between all provided subjects subsets
        * ...
        * z_value: number of standard deviations by which data point is above the mean value.
        * The considered data point is the difference between two linear regressions slopes.
        * abs_z_value: absolute z_value
        * p_value: probability of rejecting the null hypothesis that the difference in slopes is zero.
        * ...

        For each subjects subset the next columns are added to the resulting table:

        * R2_***: determination coefficient.
          A statistical measure of how well the regression line approximates the data points.
        * intercept_***: estimated value of the intercept of linear regression.
        * slope_***: estimated value of the slope of linear regression.
        * intercept_std_***: standard error of the estimate of the intercept of linear regression.
        * slope_std_***: standard error of the estimate of the slope of linear regression.
        * intercept_p_value_***: p-value for the intercept of linear regression.
        * slope_p_pvalue_***: p-value for the slope of linear regression.
        * ...

        Where *** is the name of subjects subset.

        Possible parameters of experiment:

        * None

        :param data: pdm.Data instance, which specifies information about dataset.
        :param annotations: pdm.Annotations instance, which specifies subset of CpGs.
        :param attributes: pdm.Attributes instance, which specifies information about subjects.
        :param observables_list: list of subjects subsets. Each element in list is dict,
         where ``key`` is observable name and ``value`` is possible values for this observable.
        :param method_params: parameters of experiment.
    """

    table_aggregator_linreg(
        DataType.betas,
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
    table_aggregator_variance(
        DataType.genes,
        data,
        annotations,
        attributes,
        observables_list,
        data_params=data_params,
    )
