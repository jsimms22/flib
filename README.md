# flib

ROAD MAP TO DO LIST:
[] Need to fix the array_printer for row and col major versions
[] Need to fix the array_builder for row and col major versions
[] Verify that matmul.h is working correctly despite row/col major
    >e.i. add a way to transpose from row to maj before computing
[] Remove buffered slots to return answer in original format
[] Get dgemm to work - currently returns an issue of reference_dgemm() does not recognize 
the existance of cblas_dgemm()
[] Reconsider the reason for needing a cpp version
[] Parallelize the code to work with OMP
[] Parallelize the code to work with MPI
[] Attempt to reconcile OMP with MPI
[x] Get Raspberry Pi Beowulf Cluster working
    >check that the fstab file has 777 rights, currently you have to use sudo mount -a to 
     reconcile NFS across the cluster
    >set static IPs
    >redistribute public keys for passwordless SSH from Pi to Pi for MPI
[] evaluate performance on cluster against BLAS/LAPACK IOT write paper


THINGS TO CONSIDER:
[] using mdspan vs std::array vs std::vector. Planning to just finish it using std::array as the container, but might convert to mdspan at some point
