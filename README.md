# Domain Decomposition

This repository is a complementary of a <a href="https://summerofhpc.prace-ri.eu/can-you-briefly-explain-the-domain-decomposition-method/">blog post</a> which I have written about the domain decomposition method. It contains an exercise about domain decomposition and its solution.

## Problem Definition

Write a program that gets the average of each cell and its neighbors(north, south, east, west, nort-west, north-east, south-west, south-east) and writes the result into a file in parallel.

Your solution should satisfy these:

* creating the initial matrix in parallel
* calculating the average of each cell's neighbors and writing the result into a file in parallel
* using at least 2 processes
* scaling with different number of processes

Getting average of every cell is considered as an iteration. You should do this iteration at least 2 times. Once you compare your results with [sample file](sample_output), you can play with the number of iterations to measure the performance.

You can ignore corner cases like having a 10x10 matrix and dividing it into 3 processes. You can assume the matrix is divided to equal squares among processes. In addition to that, you do not have to decompse the domain into squares to solve the problem. Any geometrical shape can be used.


The final result should look like this:
```

#Initial matrix
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0
0 0 0 1 1 0 0 0
0 0 0 1 1 0 0 0
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0


#First iteration
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0
0 0 0.111111111111111 0.222222222222222 0.222222222222222 0.111111111111111 0 0
0 0 0.222222222222222 0.444444444444444 0.444444444444444 0.222222222222222 0 0
0 0 0.222222222222222 0.444444444444444 0.444444444444444 0.222222222222222 0 0
0 0 0.111111111111111 0.222222222222222 0.222222222222222 0.111111111111111 0 0
0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0

#Second iteration
0 0 0 0 0 0 0 0
0 0.0123456790123457 0.037037037037037 0.0617283950617284 0.0617283950617284 0.037037037037037 0.0123456790123457 0
0 0.037037037037037 0.111111111111111 0.185185185185185 0.185185185185185 0.111111111111111 0.037037037037037 0
0 0.0617283950617284 0.185185185185185 0.308641975308642 0.308641975308642 0.185185185185185 0.0617283950617284 0
0 0.0617283950617284 0.185185185185185 0.308641975308642 0.308641975308642 0.185185185185185 0.0617283950617284 0
0 0.037037037037037 0.111111111111111 0.185185185185185 0.185185185185185 0.111111111111111 0.037037037037037 0
0 0.0123456790123457 0.037037037037037 0.0617283950617284 0.0617283950617284 0.037037037037037 0.0123456790123457 0
0 0 0 0 0 0 0 0
```


MPI Topics that can be useful to solve the exercise:

* Domain Decomposition
* Cartesian Topology
* Derived Datatypes
* File I/O
* File Views

## Solution

You can find the solution in C language [here](blurring_effect.c). Note that this is not the only solution. There may be dozens of different approaches. Mine is just one of them. The program runs on 4 processes. If you'd like to play with the number of processes or the domain size, you can change these macros:

```
GLOBAL_GRID_SIZE
LOCAL_GRID_SIZE
```

For instance, if you have 10 processes and the domain of a 200x200 matrix, your variables should be like:
```
#define GLOBAL_GRID_SIZE 200
#define LOCAL_GRID_SIZE 20
```
So that each process is responsible of a 20x20 matrix.

To run the program, you should have MPI library install on your system. If you have MPI installed, run below commands to execute the code:
```
git clone git@github.com:dogukanteber/domain-decomposition.git
cd domain-decomposition/
mpicc blurring_effect.c -o b && mpirun -np 4 b
```
