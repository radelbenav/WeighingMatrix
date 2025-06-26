# WeighingMatrix
A Sage package for generating integral weighing matrices. 

1) The main file is Classification.sage and the main function is ExhaustiveListIW(n,w). This function outputs the list of all integer weighing matrices of size n and weight w up to Hadamard equivalence. Some repetitions in an isomorphism class may appear. See the file for more details.

2) The file NSOKS.sage consists of a function Nsoks(n,r), which outputs all way to write n as a sum of r squares. This is faster than previous implementations.

3) The jupyter sagemath notebook Demosntration.ipynb shows how to use the package.

For both, the main reference paper is "AN ALGORITHM FOR CONSTRUCTING AND CLASSIFYING
 THE SPACE OF SMALL INTEGER WEIGHING MATRICES, https://arxiv.org/pdf/2304.09495".

