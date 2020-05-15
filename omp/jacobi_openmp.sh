#! /bin/bash
#
gcc -c -Wall -fopenmp jacobi_openmp.c
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gcc -fopenmp -o jacobi_openmp jacobi_openmp.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm jacobi_openmp.o
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
time ./jacobi_openmp > jacobi_openmp.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
time ./jacobi_openmp >> jacobi_openmp.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
time ./jacobi_openmp >> jacobi_openmp.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Run with 8 threads."
export OMP_NUM_THREADS=8
time ./jacobi_openmp >> jacobi_openmp.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
rm jacobi_openmp
#
echo "Normal end of execution."
