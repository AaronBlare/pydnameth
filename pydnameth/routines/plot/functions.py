import plotly.graph_objs as go
from pydnameth.config.experiment.types import DataType, Method
from pydnameth.routines.common import get_axis, get_legend, get_margin


def get_layout(config):
    layout = None

    if config.experiment.data in [DataType.betas, DataType.residuals_common, DataType.residuals_special,
                                  DataType.genes]:

        if config.experiment.method in [Method.scatter, Method.variance_histogram]:

            item = config.experiment.method_params['item']

            if config.experiment.data == DataType.genes:
                aux_str = ''
            else:
                if item in config.cpg_gene_dict:
                    aux = config.cpg_gene_dict[item]
                    if isinstance(aux, list):
                        aux_str = '(' + ';'.join(aux) + ')'
                    else:
                        aux_str = '(' + str(aux) + ')'
                else:
                    aux_str = '(non-genic)'

            y_title = 'Methylation level'
            if config.experiment.data in [DataType.residuals_common, DataType.residuals_special]:
                y_title = 'Residuals'

            layout = go.Layout(
                title=dict(
                    text=item + aux_str,
                    font=dict(
                        family='Arial',
                        size=33,
                    )
                ),
                autosize=True,
                margin=go.layout.Margin(
                    l=95,
                    r=10,
                    b=80,
                    t=85,
                    pad=0
                ),
                barmode='overlay',
                legend=dict(
                    font=dict(
                        family='Arial',
                        size=16,
                    ),
                    orientation="h",
                    x=0.33,
                    y=1.11,
                ),
                xaxis=get_axis(config.attributes.target.capitalize()),
                yaxis=get_axis(y_title),
            )

        elif config.experiment.method in [Method.curve]:

            x_title = config.experiment.method_params['x']
            y_title = config.experiment.method_params['y']

            layout = go.Layout(
                autosize=True,
                margin=go.layout.Margin(
                    l=95,
                    r=10,
                    b=80,
                    t=10,
                    pad=0
                ),
                barmode='overlay',
                legend=dict(
                    font=dict(
                        family='Arial',
                        size=16,
                    ),
                    orientation="h",
                    x=0.33,
                    y=1.11,
                ),
                xaxis=dict(
                    title=x_title,
                    zeroline=False,
                    showgrid=True,
                    showline=True,
                    mirror='ticks',
                    titlefont=dict(
                        family='Arial',
                        size=33,
                        color='black'
                    ),
                    showticklabels=True,
                    tickangle=0,
                    tickfont=dict(
                        family='Arial',
                        size=30,
                        color='black'
                    ),
                    exponentformat='e',
                    showexponent='all'
                ),
                yaxis=get_axis(y_title),
            )

    elif config.experiment.data == DataType.epimutations:

        if config.experiment.method == Method.scatter:

            y_title = 'SEM number'

            layout = go.Layout(
                title=dict(
                    font=dict(
                        family='Arial',
                        size=33,
                    )
                ),
                autosize=True,
                margin=go.layout.Margin(
                    l=120,
                    r=10,
                    b=80,
                    t=50,
                    pad=0
                ),
                barmode='overlay',
                legend=dict(
                    font=dict(
                        family='Arial',
                        size=16,
                    ),
                    orientation="h",
                    x=0.33,
                    y=1.11,
                ),
                xaxis=get_axis(config.attributes.target.capitalize()),
                yaxis=get_axis(y_title),
            )

        elif config.experiment.method == Method.range:

            y_title = 'SEM number'

            layout = go.Layout(
                title=dict(
                    font=dict(
                        family='Arial',
                        size=33,
                    )
                ),
                autosize=True,
                margin=go.layout.Margin(
                    l=120,
                    r=10,
                    b=70,
                    t=10,
                    pad=0
                ),
                barmode='overlay',
                showlegend=False,
                xaxis=dict(
                    title=config.attributes.target.capitalize(),
                    showgrid=True,
                    showline=True,
                    mirror='ticks',
                    titlefont=dict(
                        family='Arial',
                        size=33,
                        color='black'
                    ),
                    showticklabels=True,
                    tickangle=0,
                    tickfont=dict(
                        family='Arial',
                        size=18,
                        color='black'
                    ),
                    exponentformat='e',
                    showexponent='all'
                ),
                yaxis=get_axis(y_title),
            )

    elif config.experiment.data == DataType.entropy:

        if config.experiment.method == Method.scatter:
            y_title = 'Entropy'

            layout = go.Layout(
                title=dict(
                    font=dict(
                        family='Arial',
                        size=33,
                    )
                ),
                autosize=True,
                margin=go.layout.Margin(
                    l=120,
                    r=10,
                    b=80,
                    t=50,
                    pad=0
                ),
                barmode='overlay',
                legend=dict(
                    font=dict(
                        family='Arial',
                        size=16,
                    ),
                    orientation="h",
                    x=0.33,
                    y=1.11,
                ),
                xaxis=get_axis(config.attributes.target.capitalize()),
                yaxis=get_axis(y_title),
            )

    elif config.experiment.data == DataType.observables:

        if config.experiment.method == Method.histogram:
            layout = go.Layout(
                autosize=True,
                margin=get_margin(),
                barmode=config.experiment.method_params['barmode'],
                legend=get_legend(),
                xaxis=get_axis(config.attributes.target.capitalize()),
                yaxis=get_axis('Count'),
            )

    return layout
