"case1.R", "case2&3.R" and "function.R" are used to get seeds list that failed to generate gs boundary because of No. of events.
"case1.R" and "case2&3.R" only differs at parameter setup.


"main.R" is used to generate initial seeds firstly and calls "case1.R", "case2&3.R" and "function.R" to get final seed list secondly.
Change iter=4000 if more simultions is needed(currently set for 2000 simulations). 
For example, if 5000 simulations are needed, iter should be set to be bigger than 5000(iter=simulation times + possible failed seeds).