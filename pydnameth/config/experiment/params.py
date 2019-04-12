from pydnameth.config.experiment.types import DataType, Method, Task


def get_default_params(config):
    params = {}

    if config.experiment.type in [DataType.betas, DataType.residuals_common, DataType.residuals_special]:

        if config.experiment.task == Task.table:

            if config.experiment.method == Method.cluster:
                params = {
                    'eps': 0.1,
                    'min_samples_percentage': 1
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
                    'details': 2
                }
            elif config.experiment.method == Method.variance_histogram:
                params = {
                    'item': 'cg01620164',
                }

    elif config.experiment.type == DataType.observables:

        if config.experiment.task == Task.plot:

            if config.experiment.method == Method.histogram:
                params = {
                    'bin_size': 1.0,
                    'opacity': 0.8,
                    'barmode': 'stack',
                    'x_range': 'auto'
                }

    elif config.experiment.type == DataType.epimutations:

        if config.experiment.task == Task.plot:

            if config.experiment.method == Method.scatter:
                params = {
                    'x_range': 'auto',
                    'y_range': 'auto',
                    'y_type': 'linear'
                }

    return params
