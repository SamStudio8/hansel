Hansel
======

A graph-inspired data structure for determining likely chains of ordered symbols from breadcrumbs of evidence.
Brother to `Gretel
<https://github.com/SamStudio8/gretel>`_.

What is it?
-----------

**Hansel** is a probabilistically-weighted, graph-inspired, novel data structure.
Hansel is designed to store the number of observed occurrences of a symbol `a` appearing at some position in space or time `i`, co-occurring with another symbol `b` at another position in space or time `j`.

One may traverse along ordered positions in time or space, each time predicting the next most likely symbol of the sequence to traverse to, given the previously selected symbols in the path.
Hansel presents a user-friendly API for managing and interacting with the data stored within.

Requirements
------------
::

    pip install numpy

Install
-------
::

    pip install hanselx

Citation
--------
Paper pending...

License
-------
Hansel and Gretel are distributed under the MIT license, see LICENSE.

