"""Legend creation with duplicate-label removal."""

from __future__ import annotations


def configure_legend(ax, spec: dict | bool | None = None, *, other_axes=()):
    if spec is False:
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()
        return None
    handles, labels = ax.get_legend_handles_labels()
    for other in other_axes:
        other_handles, other_labels = other.get_legend_handles_labels()
        handles.extend(other_handles); labels.extend(other_labels)
    unique = {}
    for handle, label in zip(handles, labels):
        if label and not label.startswith("_") and label not in unique:
            unique[label] = handle
    if not unique:
        return None
    options = {} if spec is None or spec is True else dict(spec)
    return ax.legend(unique.values(), unique.keys(), **options)
