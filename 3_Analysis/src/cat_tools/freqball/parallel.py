"""Optional MPI with collective error reporting and a single output writer."""

import os


def context():
    try:
        from mpi4py import MPI
    except (ImportError, OSError, RuntimeError) as exc:
        sizes = ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMIX_SIZE", "MV2_COMM_WORLD_SIZE")
        if any(int(os.environ.get(key, "1")) > 1 for key in sizes):
            raise RuntimeError("MPI launcher detected but mpi4py could not initialize") from exc
        return None, 0, 1
    comm = MPI.COMM_WORLD
    return comm, comm.Get_rank(), comm.Get_size()


def collectively(function, comm):
    """Finish rank-local work before propagating any error to every rank."""
    if comm is None:
        return function()
    value, error = None, None
    try:
        value = function()
    except Exception as exc:  # noqa: BLE001 - synchronize failures before the next MPI collective
        error = f"rank {comm.Get_rank()}: {type(exc).__name__}: {exc}"
    errors = comm.allgather(error)
    if any(errors):
        raise RuntimeError("FreqBall parallel stage failed: " + "; ".join(x for x in errors if x))
    return value
