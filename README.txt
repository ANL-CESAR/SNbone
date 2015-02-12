Apply A mini-app version 1.0
This program emulates the inversion of A in A*x=S on a single node arch. 
This operation was identified as the bulk of the flop work and the largest memory issue for deterministic neutronics

*----------------------------------------------------------------------------------------------------------------------*
Compilation
*----------------------------------------------------------------------------------------------------------------------*
There is a Makefile in each subdirectory.
The Makefile has sufficient comments at the top to describe the options for compiling and we defer there for further instructions.
The existing makefile detects the platform type via "uname" or "hostname" identification which we suggest the user alter for their needs. 

*----------------------------------------------------------------------------------------------------------------------*
Code changes
*----------------------------------------------------------------------------------------------------------------------*
You are free do alter this code however you see fit, but note that 
1) Replacing FGMRES with other solvers (BiCGSTAB) is not advisable
2) Adding a spatial multi-grid or algebraic multi-grid step will help performance but will also alter the performance compared with the app.
   Our goal is just to understand the performance characteristics of the given subroutines with regard to architecture assumptions

*----------------------------------------------------------------------------------------------------------------------*
Program design and mini-app goals
*----------------------------------------------------------------------------------------------------------------------*
The code is broken into three steps.
1) You create an unstructured mesh but it has a regular stencil such that the bandwidth is not as severe as a conventional FE mesh
   You can theoretically provide a linear tetrahedral finite element mesh from a conventional mesh generator as an alternative. Ask for details.
2) You process this mesh with respect to the maximum number of threads you want to investigate
   The mesh will be re-ordered with respect to element and vertex such that it can be applied in n*thread independent steps.
3) You run the mini-app (fortran or c version) to look at the overall performance on a particular architecture

! --------
Depending upon the parallel machine setup, we want to vary the number of elements+angles assigned to a given process
If we build a machine with 1024 nodes where each node has 1024 cores, then we will want to deploy lots of angles per node and few elements
If we build a machine with say 16384 nodes where each node has 64 cores, then we will want to deploy lots of elements per node and few angles
So in effect, if you define an architecture, we want to understand how changing the space-angle partitioning affects the performance of the mini-app
! ------
This mini-app has three solution modes:
The AVE-1 scheme refers to the standard matrix-free matrix-vector applications where gaussian integration is used
The AVE-2 scheme refers to a slightly different approach from AVE-1 where intermediate vectors are used.
The AVE-3 scheme refers to a conventional assembled matrix in compressed sparse row storage.
! ---
To date, AVE-3 consistently outperforms AVE-1 and AVE-2 such that fully assembling the matrix and performing the solve is more efficient than
 doing the gaussian integration solve. There are clear reasons why this occurs.
 If there is a change in this behavior on a particular arch, then we really would like to hear about it.

The current application on Blue Gene technology assigns 2000-4000 vertices per process and 2-4 angles per MPI process. 
16 MPI processes are used on BlueGene/Q which has 16 physical cores per node.
In a threaded version of the code we would prefer to assign 32 times the current number of vertices per node.
This minimizes the net communication in our current algorithm and thus should yield better performance.

The suggested problems to study on a single core of BGQ under a 1 MPI per process application are
 makemesh.x    10 10 10 0
 processmesh.x 1  1
 SNaCFE.x      0  100 30 2 1
=> This will generate a mesh with 12000 elements, 2331 vertices
=> You can use this to get a base line for non-latency masking arch that is only the "best it will ever be able to do"
[SN-KERNEL] Method           Time       Est. Mem(MB)   Est. GFlops   Est. GFlops/s   GFlops/GByte
=SN-KERNEL]   AVE-1           1.00464          7.030         0.7715          0.768         112.37
=SN-KERNEL]   AVE-2           0.69231          7.030         0.7715          1.114         112.37
=SN-KERNEL]   AVE-3           0.01917          2.932         0.0007          0.036           0.24

For a 16 way threaded application on the same machine we would want to test more spatial vertices (16 times):
 makemesh.x    26 26 26 0
 processmesh.x 1  16
 SNaCFE.x      0  100 30 2 16
=> This will generate a mesh with 210912 elements, 37259 vertices
=> This will demonstrate the memory bandwidth limited performance
[SN-KERNEL] Method           Time       Est. Mem(MB)   Est. GFlops   Est. GFlops/s   GFlops/GByte
=SN-KERNEL]   AVE-1           1.38358        119.972        15.4954         11.200         132.26
=SN-KERNEL]   AVE-2           2.94123        119.972        15.4954          5.268         132.26
=SN-KERNEL]   AVE-3           0.03755         47.371         0.0135          0.358           0.29

Comparing these two, the number of iterations went up slightly (33 to 34 -> 34/33=1.03)
The total time went up by 1.37, 4.2, 1.96 where the last one shows the actual limitation of the chip as we have saturated the memory bandwidth.
We expect these to vary considerably on different arch.

For a 16 way threaded application with assignment of more angles rather than vertices we would do:
 makemesh.x    10 10 10 0
 processmesh.x 1  16
 SNaCFE.x      0  100 30 32 16
=> This will generate a mesh with 12000 elements, 2331 vertices
[SN-KERNEL] Method           Time       Est. Mem(MB)   Est. GFlops   Est. GFlops/s   GFlops/GByte
=SN-KERNEL]   AVE-1           1.16669         39.577        12.3436         10.580         319.37
=SN-KERNEL]   AVE-2           0.53785         39.577        12.3436         22.950         319.37
=SN-KERNEL]   AVE-3           0.03138         46.529         0.0111          0.354           0.24

Comparing this with the previous two, we
1) see a degradation in AVE-1 compared with the first example, but not as bad as the second case.
2) see an improvement in AVE-2 which has to do with the way the matrix-vector products are carried out.
3) see nearly identical negative benefits for AVE-3 with the second and third tests, but this is the same memory bandwidth saturation.

# -------
# What are the realistic limits?
#  We do not expect more than 2048 angles in the whole domain and really no more than 128 assigned to one node
#  We expect 1 billion vertices in the mesh at the crazy end but realistically we only need 10 million for most problems.
#  We are not interested in less than 1000 vertices per node cases as the current preconditioner breaks down.
# -------

*----------------------------------------------------------------------------------------------------------------------*
Execution
*----------------------------------------------------------------------------------------------------------------------*
# --------
# Build a mesh
# --------
#[SN-KERNEL] Version 1.0 Linear tetrahedral mesh generation...........
#[SN-KERNEL] Usage:   makemesh.x   X    Y    Z  BCmodel...............
#[SN-KERNEL] Example: makemesh.x  10   20   10  0      ...............
#[SN-KERNEL] X, Y, Z specifies the number of X-Y-Z meshes (E=X*Y*Z*12)
#[SN-KERNEL] BC 0-6  specifies how many box surfaces have vacuum b.c.s
makemesh.x 4 4 4 0 

# --------
# Process the mesh for a given number of threads
# --------
#[SN-KERNEL] Version 1.0 mesh processing utility.....................
#[SN-KERNEL] Usage:   processmesh.x  Partition Threads...............
#[SN-KERNEL] Example: processmesh.x  1         32     ...............
#[SN-KERNEL] Scheme   Use (1) METIS, (other) GREEDY..................
#[SN-KERNEL] Threads  Maximum threads to setup for in the mini-app...
# Partition refers to the decomposition used to load balance the threads
# Partition #1 refers to a METIS 4.0 approach to assign vertices and then elements to each thread
# Partition #2 refers to "GREEDY" approach where the element adjacency is used to color the grid and define thread sized pieces of work.
# --------
processmesh.x  1  32

# --------
# Run the mini-app
# --------
#[SN-KERNEL] Version 1.0 SNaCFE mini-app to study on node performance
#[SN-KERNEL] This mini-app is a test of the within-group FGMRES solver for a CFE SUPG based SN methodology
#[SN-KERNEL] Usage:   snacfe.x  Scheme Iter BackV Angles Threads
#[SN-KERNEL] Example: snacfe.x  1      100  30    32     1      
#[SN-KERNEL] Scheme          specifies which scheme to use for the study (0=all)...........
#[SN-KERNEL] Iter(ation)     specifies the maximum FGMRES iterations to allow..............
#[SN-KERNEL] Back V(ectors)  specifies the maximum back vectors to use in FGMRES...........
#[SN-KERNEL] Angles          specifies the number of angles assigned to the local process..
#[SN-KERNEL] T(hreads)       specifies the number of threads to use during the execution...
# --------
SNaCFE.x 0 100 30 8 8

*----------------------------------------------------------------------------------------------------------------------*
Example Output
*----------------------------------------------------------------------------------------------------------------------*
[SN-KERNEL] Version 1.0 SNaCFE mini-app to study on node performance....................................................
[SN-KERNEL] This mini-app is a test of the within-group FGMRES solver for a CFE SUPG based SN methodology...............
[SN-KERNEL].............................................................................................................
[SN-KERNEL] Usage:   snacfe.x  Scheme Iter BackV Angles Threads....................................................
[SN-KERNEL] Example: snacfe.x  1      100  30    32     1      ....................................................
[SN-KERNEL].............................................................................................................
[SN-KERNEL].............................................................................................................
[SN-KERNEL] Thread id     0 of     1
[SN-KERNEL] Running schemes  1: 3
[SN-KERNEL] Number of Iterations         100
[SN-KERNEL] Number of Back Vectors        30
[SN-KERNEL] Number of Angles               4
[SN-KERNEL] Number of Threads              1
[SN-KERNEL] Importing the processed mesh
[SN-KERNEL] Setting up the angle cubature
[SN-KERNEL] Stenciling the spatial NZ matrix
[SN-KERNEL] Number of Elements             2592
[SN-KERNEL] Number of Vertices              559
[SN-KERNEL] Vector Size Assembled          2236
[SN-KERNEL] Vector Size Dis-assem         41472
[SN-KERNEL] Number of Non-zeros            7291
[SN-KERNEL] Average Connections/Vertex       13
[SN-KERNEL] Allocating FGMRES memory and solution vectors
[SN-KERNEL] Building the non-zero space-angle matrices
[SN-KERNEL] Generating an answer and its associated source with size      2236
[SN-KERNEL]...Calling FGMRES solver for AVE1
[SN-KERNEL]...FGMRES returned with an error of  6.253887E-06 after     28 iterations
[SN-KERNEL] ***SUCCESS*** check passed for     AVE-1   
[SN-KERNEL]...Calling FGMRES solver for AVE2
[SN-KERNEL]...FGMRES returned with an error of  6.253880E-06 after     28 iterations
[SN-KERNEL] ***SUCCESS*** check passed for     AVE-2   
[SN-KERNEL]...Calling FGMRES solver for AVE3
[SN-KERNEL]...FGMRES returned with an error of  6.253849E-06 after     28 iterations
[SN-KERNEL] ***SUCCESS*** check passed for     AVE-3   
[SN-KERNEL] Average FGMRES memory (MB)         0.979
[SN-KERNEL] Assembly Time =     0.002383 seconds or           9.4 FGMRES iterations 
[SN-KERNEL] Method           Time       Est. Mem(MB)   Est. GFlops   Est. GFlops/s   GFlops/GByte
=SN-KERNEL]   AVE-1           0.33811          2.028         1.4279          4.223         721.07
=SN-KERNEL]   AVE-2           0.17785          2.028         1.4279          8.029         721.07
=SN-KERNEL]   AVE-3           0.00708          1.317         0.0009          0.125           0.69
[SN-KERNEL].............................................................................................................

*----------------------------------------------------------------------------------------------------------------------*
Useful stuff
*----------------------------------------------------------------------------------------------------------------------*
export KMP_AFFINITY=scatter or KMP_AFFINITY=compact are OpenMP options useful for testing node performance

