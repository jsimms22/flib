# flib

ROAD MAP TO DO LIST:
-Need to fix the array_printer for row and col major versions
-Need to fix the array_builder for row and col major versions
-Verify that matmul.h is working correctly despite row/col major
    >e.i. add a way to transpose from row to maj before computing
-Remove buffered slots to return answer in original format
-Get dgemm to work - currently returns an issue of reference_dgemm() does not recognize 
the existance of cblas_dgemm()
-Reconsider the reason for needing a cpp version
-Parallelize the code to work with OMP
-Parallelize the code to work with MPI
-Attempt to reconcile OMP with MPI
-Get Raspberry Pi Beowulf Cluster working
    >check that the fstab file has 777 rights, currently you have to use sudo mount -a to 
     reconcile NFS across the cluster
    >set static IPs
    >redistribute public keys for passwordless SSH from Pi to Pi for MPI
-evaluate performance on cluster against BLAS/LAPACK IOT write paper



01AUG20 CHANGE LOG:
1) DELETED:
  ->column major means of printing and buffing arrays in test.c (waste of effort).
  ->it will be easier to write code to print/buff inside the native language for column
    major languages.

2) NOTES:
  ->matmul works for column major situations.
  ->need to decide if I want to transpose row major arrays or allow for row major math
    using avx256 (time vs optimization issue).
  ->row_naive_dgemm works for math checking if I decide to transpose inside matmul.h or
    write optimized row major operations.
