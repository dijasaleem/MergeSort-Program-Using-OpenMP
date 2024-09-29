#include <iostream>
#include <cstdlib>
#include <ctime>
#include <mpi.h>
#include <omp.h>
#include <string>

using namespace std;

void merge(int* arr, int l, int m, int r) {
    int i, j, k;
    int s1 = m - l + 1;
    int s2 = r - m;
    int* a1 = new int[s1];
    int* a2   = new int[s2];

    for (i = 0; i < s1; i++)
        a1[i] = arr[l + i];
    for (j = 0; j < s2; j++)
        a2[j] = arr[m + 1 + j];

    i = 0;
    j = 0;
    k = l;

    while (i < s1 && j < s2) {
        if (a1[i] >= a2[j]) {
            arr[k] = a1[i];
            i++;
        }
        else {
            arr[k] = a2[j];
            j++;
        }
        k++;
    }

    while (i < s1) {
        arr[k] = a1[i];
        i++;
        k++;
    }

    while (j < s2) {
        arr[k] = a2[j];
        j++;
        k++;
    }
    delete[] a1;
    delete[] a2;
}

void mergeSort(int* arr, int left, int right) {
    if (left < right) {
        int mid = (left + right) / 2;

#pragma omp task
            {
                mergeSort(arr, left, mid);
            }
#pragma omp task
            {
                mergeSort(arr, mid + 1, right);
            }
#pragma omp taskwait       
        merge(arr, left, mid, right);
    }
}
int main(int argc, char** argv) {
    int m = 5; 
    int n = 4; //default value
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    m = size;
    srand(time(NULL));
    
    int* A = new int [m*n];
    int* row = new int[n];

    if (rank == 0) {
        // Initialize matrix randomly
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i * n + j] = rand() % 10;
            }
        }
        cout << "Original Matrix:" << endl;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                cout << A[i * n + j] << " ";
            }
            cout << endl;
        }
    }
    
    MPI_Scatter(A, n, MPI_INT, row, n, MPI_INT, 0, MPI_COMM_WORLD);
    mergeSort(row, 0, n - 1);
    cout << "I am process " << rank << endl;
    for (int i = 0; i < n; i++) {
        cout << row[i] << " ";
    }
    cout << endl;
    int sum = 0;
    for (int i = 0; i < n; i++) {
        sum += row[i];
    }
    int* sortedMatrix = NULL;
    if (rank == 0) {
        sortedMatrix = new int[m * n];
    }
    MPI_Gather(row, n, MPI_INT, sortedMatrix, n, MPI_INT, 0, MPI_COMM_WORLD);
    int totalSum = 0;
    MPI_Reduce(&sum, &totalSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Sorted matrix:" << endl;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                cout << sortedMatrix[i * n + j] << " ";
            }
            cout << endl;
        }
        cout << "Total sum: " << totalSum << endl;
    }

    //deleting dynammic memory
    delete[] A;
    delete[] row;
    if(rank == 0)
        delete[] sortedMatrix;
    MPI_Finalize();
    return 0;
}
