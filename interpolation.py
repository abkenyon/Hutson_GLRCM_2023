import interp
import numpy as np

def accadia_f(points, values, xi, x_deltas, y_deltas, query):
    """Interpolate using a subgrid nearest neighbors method (Accadia, 2003).

    This function calls a Fortran subroutine that is compiled with `f2py`. Specifically,
    at the command line run

    $ f2py -c -m interp interp.f

    which will create a python module called interpolation.

    Parameters
    ----------
    points : N x 2 array
        Coordinates of original data points.
    values : N array
        Data values at each of the coordinates.
    xi : M x 2 array
        Coordinates to analyze/interpolate to.
    x_deltas : M x 1 array
        X grid spacing of target grid.
    y_deltas : M x 1 array
        Y grid spacing of target grid.
    query : list
        Neighboring `points` for each `xi`. This should be the output of a search from a
        method as `scipy.spatial.cKDTree.query_ball_point`.

    Returns
    -------
    numpy.ndarray
    """
    n_neighbors_max = max(len(q) for q in query)
    query_arr = np.full((len(query), n_neighbors_max), -1, int)

    # Convert query to an array so it can be operated on in fortran
    for i, q in enumerate(query):
        query_arr[i, :len(q)] = np.array(q) + 1

    analysis = interp.accadia(points, values, xi, x_deltas, y_deltas, query_arr)

    return analysis

