# SNbone
This application targets the primary computational solve burden of a SN,
continuous finite element based transport equation solver

## Branch Features

### Threads vs. Task Groups

The main purpose of this branch is to decouple (a) the number of
non-overlapping work elements and (b) the number of threads.  To decouple the
two, we introduce "task groups" that are independent of the number of threads.  

In the original workflow, the mesh is preprocessed for a specified number of
threads (AS_NumThreads).  Then in SNBone, the mesh data are distributed across
the available threads (TasksPerThread, AS_ThreadWiseWork, etc.).  The conditions
are as follows.  

   ```
   AS_NumThreads >= NumThreads
   TasksPerThread = AS_NumThreads / NumThreads
   ```

With the introduction of task groups, the mesh is still be preprocessed for a
specified number of threads.  However, in SNBone, the user will specify a
number of task groups (NumTaskGroups) and then the mesh data are distributed
across task groups with the following conditions:

   ```
   AS_NumThreads >= NumTaskGroups
   TasksPerGroup = As_NumThreads / NumTaskGroups
   ```

Each task group has specified start- and end-points in the data decomposition
schemes, as threads do in the original implementation.  The parallelism
is be expressed slightly differently in SNBone:

```fortran
!$omp parallel

   !$omp do reduction(+:foo)
   do i=1,NumTaskGroups
      !    Process task groups in parallel, using 
      !    omp reductions when necessary.
   enddo
   !$omp end do
   ! ...implicit barrier...

   !$omp single
   !   When necessary, do a serial operation.
   !$omp end single
   ! ...implicit barrier...

   !$omp do reduction(+:bar)
   do i=1,NumTaskGroups
   !   ..etc...
   !$omp end do
   ! ...implicit barrier...

!$omp end parallel
```

For models that expose the number of threads to the user (OpenMP, CUDA), the
user can still specify the number of threads.  With OpenMP on CPU,  the new
version should hopefully match the performance of the original version when
NumThreads = NumTaskGroups (though this is unconfirmed as of writing this
README).  

For models that do not expose the number of threads to the user (OpenACC and
task-based programming models), the new task-group semantics should be
advantageous, as the number of threads is arbitrarily independent of the number
of task-groups.  However, performance will likely rely on tuning the number of
and size of task-groups for a given programming model and hardware.  
