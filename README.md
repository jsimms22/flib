# flib

##ROAD MAP TO DO LIST:
- [x] Need to fix the array_printer for row
- [x] Need to fix the array_builder for row
- [x] Verify that matmul.h is working correctly 
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


##SAMPLE OUTPUT
Beginning <double> std::array tests

C matrix solution from function naive_row_matmul:
Time taken by function: 1977700 microseconds
C solution from function avx256_row_matmul:
Time taken by function: 828900 microseconds

Beginning <double> heap std::array pointer tests

C matrix solution from function naive_row_matmul:
Time taken by function: 2324700 microseconds
C solution from function avx256_row_matmul:
Time taken by function: 992200 microseconds

Beginning <float> std::array tests

C matrix solution from function naive_row_matmul:
Time taken by function: 2887500 microseconds
C solution from function avx256_row_matmul:
Time taken by function: 557000 microseconds