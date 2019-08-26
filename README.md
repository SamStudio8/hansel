<div align="center">
<p align="center">
    <img src="hansel-logo.png?raw=true?" alt="gretel-logo" width="200">
</p>
<h1 align="center">Hansel</h1>
<h3 align="center">A graph-inspired data structure for determining likely chains of ordered symbols from breadcrumbs of evidence. Brother of <a href="https://github.com/SamStudio8/gretel">Gretel</a>.
</h3>
<p align="center">
<a href="https://github.com/samstudio8/hansel/blob/master/LICENSE"><img src="https://img.shields.io/badge/license-MIT-orange.svg" alt="License"></a>
<a href="https://bioconda.github.io/recipes/hanselx/README.html"><img src="https://anaconda.org/bioconda/hanselx/badges/downloads.svg" alt="bioconda"></a>
</p>
</div>

What is it?
-----------

**Hansel** is a probabilistically-weighted, graph-inspired, novel data structure.
Hansel is designed to store the number of observed occurrences of a symbol `a` appearing at some position in space or time `i`, co-occurring with another symbol `b` at another position in space or time `j`.

One may traverse along ordered positions in time or space, each time predicting the next most likely symbol of the sequence to traverse to, given the previously selected symbols in the path.
Hansel presents a user-friendly API for managing and interacting with the data stored within.

Requirements
------------

    pip install numpy

Install
-------

    pip install hanselx

Citation (Pre-Print)
--------
```
@article {Nicholls223404,
	author = {Nicholls, Samuel M. and Aubrey, Wayne and Edwards, Arwyn and de Grave, Kurt and Huws, Sharon and Schietgat, Leander and Soares, Andr{\'e} and Creevey, Christopher J. and Clare, Amanda},
	title = {Computational haplotype recovery and long-read validation identifies novel isoforms of industrially relevant enzymes from natural microbial communities},
	elocation-id = {223404},
	year = {2018},
	doi = {10.1101/223404},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2018/01/13/223404},
	eprint = {https://www.biorxiv.org/content/early/2018/01/13/223404.full.pdf},
	journal = {bioRxiv}
}
```

License
-------
Hansel and Gretel are distributed under the MIT license, see LICENSE.

