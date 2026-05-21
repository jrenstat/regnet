# regnet 1.1.0

* Added `foldid` to `cv.regnet()` for user-controlled cross-validation folds.
* Added an adjacency construction option for choosing thresholded or full adjacency matrices.
* Added an S3 method for fitting `regnet()` from a `cv.regnet` object and explicitly supplied lambda values.
* Updated package metadata and README source.

# regnet 1.0.2 [2025-2-9]

* Bug Fix: Fixed a bug in the C++ function `QRWMR`.
* Maintenance: Updated Maintainer email.


# regnet 1.0.0 [2022-8-12]

* Added the robust network regularization for the continuous response.
* A generic function plot() is added for plotting the network structures among the identified genetic variants.

# regnet 0.4.0 [2019-6-7]

Based on users' feedback, we have
* Added more checking steps for data format, which help users make sure their data are in the correct format.
* Provided more information in the documentation for troubleshooting .

# regnet 0.3.0 [2018-5-21]

* Two new, easy to use, integrated interfaces: cv.regnet() and regnet().
* New methods for continuous and survival responses.
* The new "clv" argument allows the presence of clinical variables that are not subject to penalty in the X matrix.

# regnet 0.2.0 [2017-10-14]

* Provide c++ implementation for coordinate descent algorithms, which significantly increased the speed of cross-validation functions in this package.

# regnet 0.1.1 [2017-05-23]




