This program processes a given linear tetrahedral finite element mesh for use with the mini-app
The basic idea is to reorder the elements into thread safe batches.
The coloring of the domain allows the assembly, disassembly, and matrix-vector operations to be carried out in a threaded manner.
Note that the matrix assembly routine in the mini-app is also threaded and follows the matrix assembly coloring
