# gated group sequential design (gGSD)
This repository contains the R programs of the gated group sequential design (gGSD) simulation. 
The proposed gated group sequential design controls the familywise error rate and allows multiple interim analyses to enable early stopping for efficacy or futility. It has the potential to save drug development cost and more quickly fulfill unmet medical needs.


#####################################################################################################
function.R includes log-rank test statistic calculation functions and alpha reallocation functions for gs, ad and ht.
iteration.R includes 2 functions: sf_iteration and fs_iteration. Parameters called for 3 methods in sf(start subgroup testing in ht) and
fs(start full group testing in ht) are different.

case1_sf.R: case 1 with hierarchical test starting with subgroup test.
case1_fs.R: case 1 with hierarchical test starting with full group test(only differs with case1_sf.R in iteration function).
case23_sf.R: case 2&3 with hierarchical test starting with subgroup test.
case23_fs.R: case 2&3 with hierarchical test starting with full group test(only differs with case23_sf.R in iteration function).
These 4 programs all call "function.R" and "iteration.R".
