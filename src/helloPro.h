//
// Created by Administrator on 2021/9/5.
//
#include "mpi.h"
#ifndef MYMPI_DEMO_HELLOPRO_H
#define MYMPI_DEMO_HELLOPRO_H
void helloPrint(int argc, char** argv);
void barrierDemo(int argc, char** argv);
void bcastDemo(int argc, char** argv);
void scatterDemo(int argc, char** argv);
void gatherDemo(int argc, char** argv);
void reduceDemo(int argc, char** argv);
void officalReduceDemo(int argc, char **argv);
void alltoallDemo(int argc, char **argv);
void nonBlockingDemo(int argc, char **argv);
void diyDemo(int argc, char **argv);
void openMpDemo(int argc, char **argv);
void quickSort(int *num,int low,int high);//进行分区
int Partition(int *num,int low,int high);//返回分离点
void read_num(long long int *num_point,int my_rank,MPI_Comm comm);
void compute_pi(long long int num_point,long long int* num_in_cycle,long long int* local_num_point,int comm_sz,long long int *total_num_in_cycle,MPI_Comm comm,int my_rank);
void PiDemo(int argc, char **argv);
#endif //MYMPI_DEMO_HELLOPRO_H
