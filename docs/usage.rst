Usage
=====

Adding Observations
-------------------

Simple
~~~~~~
To construct and add observations to a **Hansel** data structure: ::

    from hansel import Hansel
    import numpy as np

    symbols = ['A', 'C', 'G', 'T', '_']
    unsymbols = ['_']
    positions = [1, 3, 5, 7, 9]
    L = 10 # Defines the 'lookback'

    # Get some memory and pass it to the Hansel constructor
    a = np.ndarray((
        len(symbols), len(symbols), len(positions)+2, len(positions)+2)
    )
    hansel = Hansel(a, symbols, unsymbols, L=L)

    # Add some observations
    hansel.add_observation(a, b, i, j)
    # ...


Not so Simple
~~~~~~~~~~~~~
For very large data sets, or complicated parallel high-throughput methodologies,
you may need to bypass the API for adding observations by reserving and populating
the relevant memory yourself: ::

    import ctypes
    import numpy as np
    from hansel import Hansel

    symbols = ['A', 'C', 'G', 'T', '_']
    unsymbols = ['_']
    positions = [1, 3, 5, 7, 9]
    L = 10 # Defines the 'lookback'

    # Get some memory
    h = np.frombuffer(Array(ctypes.c_float,
            (len(symbols)**2) * ((len(positions)+2)**2),
            lock=False),
    dtype=ctypes.c_float)

    # Shape the memory into a numpy array of the desired size
    h = hansel.reshape(len(symbols), len(symbols), len(positions)+2, len(positions)+2)
    h.fill(0.0)

    # Add some observations
    def __symbol_num(symbol):
        symbols_d = {symbol: i for i, symbol in enumerate(symbols)}
        return symbols_d[symbol]
    h[__symbol_num(a), __symbol_num(b), i, j] += 1
    # ...

    # Feed the prefilled array to the Hansel constructor
    hansel = Hansel(h, symbols, unsymbols, L=L)

Get Observation Counts
----------------------

To find the number of times symbol `a` at position `i` has been seen with symbol
`b` at position `j`: ::

    hansel.get_observation(a, b, i, j)

Summarise Counts at Position
----------------------------

To fetch a dictionary of the raw counts for each symbol at a given position: ::

    hansel.get_counts_at(at_position)

Marginal Distribution
---------------------

To find the `log10` probability of a particular symbol
appearing at a given position: ::

    hansel.get_marginal_of_at(interesting_symbol, at_position)

Conditional Distribution
------------------------
To find the `log10` probability of `a` at `i` appearing with `b` and `j`: ::

    hansel.get_conditional_of_at(a, b, i, j)

Get Spanning Support
--------------------
Get the number of times a symbol or state `b` appears at position `j`, on pieces of evidence that also covered
space or time point `i`: ::

    hansel.get_spanning_support(b, i, j)

Get Edge Weights
----------------
Given a sequence of symbols selected during traversal thus far, find the `log10` probabilities of traversing to the available symbols at your next position `j`: ::

    hansel.get_edge_weights_at(j, current_path)

Reweight Observations
---------------------
Reduce the element that support the observation of `a` at `i` and `b` at `j` co-occurring on the same piece of evidence together: ::

    hansel.reweight_observation(a, b, i, j, ratio)

The original observation count is multiplied by the ratio, the result is then subtracted from the current value.
It is recommended that `ratio` not be too large without good confidence.
Aggressive reweighting can lead to spending (removing) the evidence in the Hansel matrix before your algorithm has had time to explore the paths properly.
