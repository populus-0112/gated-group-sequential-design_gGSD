# gated group sequential design (gGSD)
This repository contains the R programs of the gated group sequential design (gGSD) simulation. 
The proposed gated group sequential design controls the familywise error rate and allows multiple interim analyses to enable early stopping for efficacy or futility. It has the potential to save drug development cost and more quickly fulfill unmet medical needs.


###########################################################################################
function.R includes log-rank test statistic calculation functions and alpha reallocation functions for GSD, AD and gGSD.


case1_sf.R: simulation for setting 1.
case23_sf.R: simulation for setting 2 and setting 3.
These 2 programs all call "function.R" and "iteration.R".
