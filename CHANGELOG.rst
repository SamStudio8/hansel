History
=======

0.0.7
-----
* Documentation

0.0.6
-----
* Switch to dictionaries and sets generated at construction time to speed
  up the lookup of symbols

0.0.5
-----
* Introduce `unsymbols` to cover cases where you would like to count observations
  from some `a` to an "invalid" `b` during marginal calculation, but still prevent
  the actual selection of the invalid `b` as a transition choice
* Member list `unsymbols` keeps track of symbols who should have no
  weighting when counting observations or calculating marginals.
* Allow a user to construct Hansel with an argument for the lookback order.
* Alter `get_spanning_support` and `get_counts_at` to not count observations
  that originate from an ignored unsymbol and transition to a real symbol.
  This alteration makes things work even if you've been a silly and filled
  Hansel with more crummy data than usual, hooray.
* Dramatically improve performance by using the correct ndarray indexing
  syntax [A, B, i, j] vs. [A][B][i][j]

0.0.4
-----
* Add some documentation.
* Rename `get_marginal_at` to `get_counts_at`. As the function returns raw
  counts, not a marginal distribution, this is a less misleading name.
* Don't return the unused `curr_branches_tot` value from `get_edge_weights_at`.
* Remove `select_next_edge_at`, we need not concern ourselves with problem
  specific end-user behaviour. We just provide an API to the pseudo-graph.

0.0.3
-----
* Add `observations` property for those who may find `crumbs` confusing or odd.
* Remove domain specific language ("SNP", "mallele") in favour of "symbol".
* Require symbol list on constuction, prevent empty list with casting/template.
* Ensure to catch an in-progress __new__ in __array_finalize__

0.0.2
-----
* Abstract BAM specific loading to `gretel`.
* Rename `reads` attribute to `slices` (of bread)
* Add `sources` property for those who may find `slices` confusing or bizarre.

0.0.1
-----
* Import repository from `claw`.
