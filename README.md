This C++ program performs matrix operations using MPI and OpenMP, including:

1. **Matrix Initialization**: Random `m x n` matrix.
2. **Process Distribution**: Each row is assigned to an MPI process.
3. **Merge Sort**: Each row is sorted in descending order using OpenMP for parallelism.
4. **Row Summation**: Each process computes the sum of its row.
5. **Result Gathering**: Process 0 gathers sorted rows and computes the total matrix sum.
