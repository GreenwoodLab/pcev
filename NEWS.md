# pcev 2.2.2

* Fixed to remove warnings due to conditions of length greater than one and recycling arrays in vector-array arithmetic

* Updated CITATION file

# pcev 2.2.1

* Added adaptive selection of blocks; see documentation for `computePCEV`.

# pcev 2.1.1

* Added `estimation = "singular"`; this option uses reduced-rank SVD for computing the component and the largest roo.

* Added methods `roysPval.PcevSingular` and `permutePval.PcevSingular`

* Expanded BLK data with multiple cell types.

# pcev 1.1.2

* Fixed the Johnstone approximation for Roy's largest root test.

* Removed VIMPblock, as it is meaningless.

* Added the possibility of computing multiple components (only for `estimation = "all"`).

* Changed how default values for `computePCEV` are handled internally.

* Throw a meaningful error when there is missing data.

# pcev 1.1.1

* Added vignette and examples.

* Changed the exact test when using shrinkage estimator.

# pcev 1.0.0

* First release.
