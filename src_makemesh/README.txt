This program creates a mesh usable for the CFE based SN mini-app.
Note that you must process this mesh before using with the mini-app.

500 FORMAT('[SN-KERNEL] Version 1.0 Linear tetrahedral mesh generation...........',51('.'))
502 FORMAT('[SN-KERNEL] Usage:   makemesh.x   X    Y    Z  BCmodel...............',51('.'))
503 FORMAT('[SN-KERNEL] Example: makemesh.x  10   20   10  0      ...............',51('.'))
509 FORMAT('[SN-KERNEL] X, Y, Z specifies the number of X-Y-Z meshes (E=X*Y*Z*12)',51('.'))
508 FORMAT('[SN-KERNEL] BC 0-6  specifies how many box surfaces have vacuum b.c.s',51('.'))


At this point, the vacuum b.c.s are not supported in the mini-app.

