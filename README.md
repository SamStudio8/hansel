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

Citation
--------
```
@article{10.1093/bioinformatics/btaa977,
    author = {Nicholls, Samuel M and Aubrey, Wayne and De Grave, Kurt and Schietgat, Leander and Creevey, Christopher J and Clare, Amanda},
    title = "{On the complexity of haplotyping a microbial community}",
    journal = {Bioinformatics},
    volume = {37},
    number = {10},
    pages = {1360-1366},
    year = {2021},
    month = {01},
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btaa977},
    url = {https://doi.org/10.1093/bioinformatics/btaa977},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/37/10/1360/38663805/btaa977.pdf},
}
```
[Read more on Twitter](https://twitter.com/samstudio8/status/1329406136592834564)


License
-------
Hansel and Gretel are distributed under the MIT license, see LICENSE.

