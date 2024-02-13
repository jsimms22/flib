# flib

##ROAD MAP TO DO LIST:
- [] Need to fix the array_printer for row
- [] Need to fix the array_builder for row
- [] Verify that matmul.h is working correctly 
- [] Remove buffered slots to return answer in original format
- [] Get dgemm to work - currently returns an issue of reference_dgemm() does not recognize the existance of cblas_dgemm()
- [] Parallelize the code to work with OMP
- [] Parallelize the code to work with MPI
- [] Attempt to reconcile OMP with MPI
- [x] Get Raspberry Pi Beowulf Cluster working
    >[x] check that the fstab file has 777 rights, currently you have to use sudo mount -a to 
     reconcile NFS across the cluster
    >[x]set static IPs
    >[x]redistribute public keys for passwordless SSH from Pi to Pi for MPI
- [] evaluate performance on cluster against BLAS/LAPACK


##THINGS TO CONSIDER:
- [] using mdspan vs std::array vs std::vector. Planning to just finish it using std::array as the container, but might convert to mdspan at some point
    > std::array or std::vector to c-style array is trivial, reverse path requires an extra step

##REFERENCES
- https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html
- https://netlib.org/lapack/explore-html/dd/d09/group__gemm_ga1e899f8453bcbfde78e91a86a2dab984.html#ga1e899f8453bcbfde78e91a86a2dab984