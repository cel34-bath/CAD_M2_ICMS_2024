This is the GitHub page for Corin Lee's talk at ICMS 2024 "Cylindrical Algebraic Decomposition in Macaulay2".

In it you will find the slides from the talk as well as some code for Macaulay2:

* CADecomposition.m2 - the package file for CADecomposition, allowing the user to perform an open CAD of R^n that is sign-invariant with respect to the input polynomials.
* RealRootsNew.m2 - an interim fix for the RealRoots package, upon which CADecomposition relies. The old RealRoots needs to be unloaded and this needs to be loaded instead.
* CADecompositionExamples.m2 - this contains some examples and things to play around with, along with commands for how to correctly load the commands.

Macaulay2 can be found at https://macaulay2.com/ with a web-based version available at https://www.unimelb-macaulay2.cloud.edu.au/.

**To run this in the Web version**
* Open Web version, hit Start.
* Open all three of the above files.
* Navigate to CADecompositionExamples.m2 and run the first few lines to install CADecomposition, uninstall the old RealRoots and install RealRootsNew.

Any questions, comments or thoughts, feel free to email me at cel34@bath.ac.uk
