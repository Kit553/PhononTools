"""Categorical-axis preparation."""

from __future__ import annotations

import pandas as pd

from catplot.data.columns import select_column


def categories(data, name, order=None):
    series = select_column(data, name)
    values = list(order) if order else list(pd.unique(series.dropna()))
    return series, values, {value: index for index, value in enumerate(values)}

