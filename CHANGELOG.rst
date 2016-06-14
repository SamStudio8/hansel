History
=======

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
