from pydnameth.config.experiment.types import DataType, Method, Task


def get_default_method_params(config):
    params = {}

    if config.experiment.data in [DataType.betas, DataType.betas_adj, DataType.residuals_common,
                                  DataType.residuals_special, DataType.genes]:

        if config.experiment.task == Task.table:

            if config.experiment.method == Method.cluster:
                params = {
                    'eps': 0.1,
                    'min_samples_percentage': 1
                }

            elif config.experiment.method == Method.variance:
                params = {
                    'semi_window': 2,
                    'box_b': 'left',
                    'box_t': 'right'
                }

            elif config.experiment.method == Method.polygon:
                params = {
                    'method': Method.linreg
                }

        elif config.experiment.task == Task.clock:

            if config.experiment.method == Method.linreg:
                params = {
                    'type': 'all',
                    'part': 0.25,
                    'size': 100,
                    'runs': 100,
                }

        elif config.experiment.task == Task.plot:

            if config.experiment.method == Method.scatter:
                params = {
                    'item': 'cg01620164',
                    'x_range': 'auto',
                    'y_range': 'auto',
                    'line': 'yes',
                    'add': 'none',
                    'semi_window': 'none'
                }
            elif config.experiment.method == Method.variance_histogram:
                params = {
                    'item': 'cg01620164',
                }
            elif config.experiment.method == Method.curve:
                params = {
                    'x': 'x',
                    'y': 'y',
                    'number_of_points': 100
                }

    elif config.experiment.data == DataType.observables:

        if config.experiment.task == Task.plot:

            if config.experiment.method == Method.histogram:
                params = {
                    'bin_size': 1.0,
                    'opacity': 0.8,
                    'barmode': 'stack',
                    'x_range': 'auto'
                }

    elif config.experiment.data == DataType.epimutations:

        if config.experiment.task == Task.plot:

            if config.experiment.method == Method.scatter:
                params = {
                    'x_range': 'auto',
                    'y_range': 'auto',
                    'y_type': 'linear',
                }

            if config.experiment.method == Method.range:
                params = {
                    'borders': [0, 20, 40, 60, 80, 100, 120],
                    'x_range': 'auto',
                    'y_range': 'auto',
                    'y_type': 'linear'
                }

    elif config.experiment.data == DataType.entropy:

        if config.experiment.task == Task.plot:

            if config.experiment.method == Method.scatter:
                params = {
                    'x_range': 'auto',
                    'y_range': 'auto',
                }

    return params
