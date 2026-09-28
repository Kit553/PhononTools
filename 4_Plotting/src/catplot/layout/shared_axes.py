"""Remove redundant labels from shared regular grids."""

from __future__ import annotations


def hide_redundant_labels(layout, *, share_x=False, share_y=False):
    for index, ax in enumerate(layout.axes.values()):
        row, column = divmod(index, layout.columns)
        if share_x and row < layout.rows - 1:
            ax.tick_params(labelbottom=False); ax.set_xlabel("")
        if share_y and column > 0:
            ax.tick_params(labelleft=False); ax.set_ylabel("")

