### Efficient Parallel Prefix/Scan in C++ using OpenMP

This code was written as part of the requirement for CSE 570 Introduction to Parallel and Distributed Programming at State University of New York at Buffalo.

The code in this repo implements the shared memory version of the popular and seminal algorithm used across parallel computing called Parallel Prefix/Scan. This code was tested and benchmarked on UB's HPC center, [Center for Computational Research.](http://www.buffalo.edu/ccr.html)

**Running the code** 

- Compiling

  - g++ -fopenmp --std==c++17 a0.cpp -O3 -o a0

- Execution

  - ./a0 n

    where n is the number of cores/processing elements

**References**

https://www.cs.cmu.edu/~guyb/papers/Ble93.pdf