//
// Created by Administrator on 2021/9/5.
//
#include "helloPro.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include "mpi.h"
#include <omp.h>
#include<math.h>
#include<time.h>

void helloPrint(int argc, char **argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int message[2];
    int dest, src;
    int tag = 0;
    MPI_Status status;
    if (size == 1) {
        printf("This example requires more than one process to execute\n");
        MPI_Finalize();
        exit(0);
    }
    if (rank != 0) {
        message[0] = rank;
        message[1] = size;
        dest = 0;
        MPI_Send(message, 2, MPI_INT, dest, tag, MPI_COMM_WORLD);
    } else {
        for (src = 1; src < size; ++src) {
            MPI_Recv(message, 2, MPI_INT, src, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            printf("Hello from process %d of %d\n", message[0], message[1]);

        }
    }
    MPI_Finalize();
}

void barrierDemo(int argc, char **argv) {
    int rank, size, len;
    MPI_Init(&argc, &argv);
    char name[MPI_MAX_PROCESSOR_NAME];

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(name, &len);

    printf("Hello,world!Process %d of %d on %s\n", rank, size, name);
    MPI_Finalize();
}

void bcastDemo(int argc, char **argv) {
    int rank, size, len;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(name, &len);
    int A[5];
//    Initialize array
    for (int i = 0; i < 5; ++i) {
        A[i] = 0;
    }
    int root = 0;
    if (rank == root) {
        A[0] = 3;
        A[1] = 4;
        A[2] = 5;
        A[3] = 2;
        A[4] = 1;
    }
    MPI_Bcast(A, 5, MPI_INT, root, MPI_COMM_WORLD);
    printf("Rank %d A[0]=%d A[1]=%d A[2]=%d A[3]=%d A[4]=%d on %s\n",
           rank, A[0], A[1], A[2], A[3], A[4], name);
    MPI_Finalize();
}

void scatterDemo(int argc, char **argv) {
    int rank, size, len;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(name, &len);
    if (size != 5) {
        printf("Example is designed for 5 processes\n");
        MPI_Finalize();
        exit(0);
    }
    int A[5], B[5];
//    initialize array
    for (int i = 0; i < 5; ++i) {
        A[i] = 0;
        B[i] = 0;
    }
//    B[0] = 0, B[1] = 0;
    int root = 0;
    if (rank == root) {
        A[0] = 3;
        A[1] = 4;
        A[2] = 5;
        A[3] = 2;
        A[4] = 1;
    }
    MPI_Scatter(A, 1, MPI_INT, B, 1, MPI_INT, root, MPI_COMM_WORLD);
    printf("Rank %d B[0] = %d B[1] = %d B[2] = %d B[3] = %d B[4] = %d on %s\n", rank, B[0], B[1], B[2], B[3], B[4],
           name);
    MPI_Finalize();
}

void gatherDemo(int argc, char **argv) {
    int rank, size, len;
    char name[MPI_MAX_PROCESSOR_NAME];
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Get_processor_name(name, &len);

        if (size != 5) {
            printf("Example is designed for 5 processes\n");
            MPI_Finalize();
            exit(0);
        }
        int A[5], B[5];

    //    initialize array
        for (int i = 0; i < 5; ++i) {
            A[i] = 0;
            B[i] = 0;
        }
        A[0] = rank;

        int root = 0;
        MPI_Gather(A, 1, MPI_INT, B, 1, MPI_INT, root, MPI_COMM_WORLD);
        printf("Rank %d B[0] = %d B[1] = %d B[2] = %d B[3] = %d B[4] = %d on %s\n", rank, B[0], B[1], B[2], B[3], B[4],
               name);
        MPI_Finalize();
    }

void reduceDemo(int argc, char **argv) {
    int rank, size, len;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(name, &len);
    int root = 0;

    if (size != 5) {
        printf("Example is designed for 5 processes\n");
        MPI_Finalize();
        exit(0);
    }
//    向量的大小为100
//    int local_vector_size = 100;
//    int n = local_vector_size;
//
//    double *b;
//    double *a;
//    a = (double *) malloc(local_vector_size * sizeof(double));
//    b = (double *) malloc(local_vector_size * sizeof(double));
//// 初始化向量的每一个元素的值
//    for (int i = 0; i < local_vector_size; ++i) {
//        a[i] = 3.14 * rank;
//        b[i] = 6.67 * rank;
//    }
//    double partial_sum = 0.0;
//// partial_sum结果是两个向量的点乘积
//    for (int i = 0; i < local_vector_size; i++) {
//        partial_sum += a[i] * b[i];
//    }
//    double sum = 0;
//    printf("Rank is %d and partial_sum = %g\n", rank, partial_sum);
//    MPI_Reduce(&partial_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);

    struct Type {
        double value;
        int location;
    } Type;
    struct Type res;
    res.location = 1;
    res.value = 0.0;
    double num = rank * 6.67;
//    printf("Rank is %d and num = %g\n",rank, num);

    MPI_Reduce(&num, &res, 2, MPI_DOUBLE_INT, MPI_MAXLOC, root, MPI_COMM_WORLD);
    if (rank == root) {
//        printf("The dot product is %g\n",sum);
//        printf("the max value is %g and location is %d\n", res.value, res.location);
    }

    printf("Rank is %d and the max value is %g and location is %d\n",rank, res.value, res.location);
//    free(a);
//    free(b);
    MPI_Finalize();
}

void officalReduceDemo(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int LEN = 1000;
    float val[LEN];        /* local array of values */
    int count=500;             /* local number of values */
    int myrank, minrank, minindex;
    float minval;
    int root = 0;

    struct {
        float value;
        int   index;
    } in, out;


/* local minloc */
    in.value = val[0];
    in.index = 0;
    for (int i=1; i < count; i++)
        if (in.value > val[i]) {
            in.value = val[i];
            in.index = i;
        }


/* global minloc */
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    in.index = myrank*LEN + in.index;
    MPI_Reduce( &in, &out, 1, MPI_FLOAT_INT, MPI_MINLOC, root, MPI_COMM_WORLD );
    /* At this point, the answer resides on process root
     */
    if (myrank == root) {
        /* read answer out
         */
        minval = out.value;
        minrank = out.index / LEN;
        minindex = out.index % LEN;
        printf("minval is %g minrank = %d minindex = %d\n", minval, minrank, minindex);
    }

    MPI_Finalize();
}

void alltoallDemo(int argc, char **argv) {
    int rank, size, len;
    char name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(name, &len);

    if (size != 5) {
        printf("Example is designed for 5 processes\n");
        MPI_Finalize();
        exit(0);
    }
    int A[5],B[5];
    for (int i = 0; i < 5; ++i) {
        A[i]=i+1+5*rank;
    }

    MPI_Alltoall(A,1,MPI_INT,B,1,MPI_INT,MPI_COMM_WORLD);
    printf("Rank: %d B: %d %d %d %d %d\n",rank,B[0],B[1],B[2],B[3], B[4]);

}

void nonBlockingDemo(int argc, char **argv) {
    int a,b;
    int size, rank;
    int tag = 0;
    MPI_Status status;
    MPI_Request send_request,recv_request;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (size != 2) {
        printf("Example is designed for 2 processes\n");
        MPI_Finalize();
        exit(0);
    }

    if (rank == 0) {
        a = 314159;
        MPI_Isend(&a,1,MPI_INT,1,tag,MPI_COMM_WORLD,&send_request);
        MPI_Irecv(&b,1,MPI_INT,1,tag,MPI_COMM_WORLD,&recv_request);

        MPI_Wait(&send_request,&status);
        MPI_Wait(&recv_request,&status);

        printf("Process %d received value %d\n",rank,b);
    } else {
        a = 677;
        MPI_Isend(&a,1,MPI_INT,0,tag,MPI_COMM_WORLD,&send_request);
        MPI_Irecv(&b,1,MPI_INT,0,tag,MPI_COMM_WORLD,&recv_request);
        MPI_Wait(&send_request,&status);
        MPI_Wait(&recv_request,&status);
        printf("Process %d received value %d\n",rank,b);
    }
    MPI_Finalize();
}

void diyDemo(int argc, char **argv) {
    typedef struct {
        int max_iter;
        double t0;
        double tf;
        double xmin;
    } Pars;

    MPI_Init(&argc,&argv);
    int rank;
    int root = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    Pars pars;
    if (rank == root) {
        pars.max_iter=10;
        pars.t0=0.0;
        pars.tf=1.0;
        pars.xmin=-5.0;
    }

    int nitems=4;

    MPI_Datatype types[nitems];

    MPI_Datatype mpi_par;
    MPI_Aint offsets[nitems];
    int blocklengths[nitems];

    types[0]=MPI_INT;offsets[0]=offsetof(Pars,max_iter);blocklengths[0]=1;
    types[1]=MPI_DOUBLE;offsets[1]=offsetof(Pars,t0);blocklengths[1]=1;
    types[2]=MPI_DOUBLE;offsets[2]=offsetof(Pars,tf);blocklengths[2]=1;
    types[3]=MPI_DOUBLE;offsets[3]=offsetof(Pars,xmin);blocklengths[3]=1;

    MPI_Type_create_struct(nitems,blocklengths,offsets,types,&mpi_par);

    MPI_Type_commit(&mpi_par);


    MPI_Bcast(&pars,1,mpi_par,root,MPI_COMM_WORLD);

    printf("Hello from rank %d; my max_iter value is %d\n",rank,pars.max_iter);

    MPI_Finalize();


}

void openMpDemo(int argc, char **argv){
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        printf("%d\n", id);
        if (id == 3) {
            printf("我是三号线程、、、");
        }
    }
}

// qwe.cpp : 定义控制台应用程序的入口点。
//





void quickSort(int *num,int low,int high)
{
    if(low<high)
    {
        int split=Partition(num,low,high);
#pragma omp parallel sections
        {
#pragma omp section
            quickSort(num,low,split-1);
#pragma omp section
            quickSort(num,split+1,high);
        }

    }

}

int Partition(int *num,int low,int high)
{
    int temp=num[low];
    while(low<high)
    {
        while(low<high&&num[high]>=temp)high--;
        num[low]=num[high];
        while(low<high&&num[low]<=temp)low++;
        num[high]=num[low];
    }
    num[low]=temp;
    return low;
}

void PiDemo(int argc, char **argv){
    long long int num_in_cycle,num_point,total_num_in_cycle,local_num_point;
    int my_rank,comm_sz;
    MPI_Comm comm;
    MPI_Init(argc,argv);
    comm=MPI_COMM_WORLD;
    MPI_Comm_size(comm,&comm_sz);//得到进程总数
    MPI_Comm_rank(comm,&my_rank);//得到进程编号
    read_num(&num_point,my_rank,comm);//读取输入数据
    compute_pi(num_point,&num_in_cycle,&local_num_point,comm_sz,&total_num_in_cycle,comm,my_rank);
    MPI_Finalize();
}

void read_num(long long int* num_point,int my_rank,MPI_Comm comm){
    if(my_rank==0){
        printf("please input num in sqaure \n");
        scanf("%lld",num_point);
    }
    MPI_Bcast(num_point,1,MPI_LONG_LONG,0,comm);

}
void compute_pi(long long int num_point,long long int* num_in_cycle,long long int* local_num_point,int comm_sz,long long int *total_num_in_cycle,MPI_Comm comm,int my_rank){
    *num_in_cycle=0;
    *local_num_point=num_point/comm_sz;
    double x,y,distance_squared;
    srand(time(NULL));
    for(long long int i=0;i< *local_num_point;i++){
        x=(double)rand()/(double)RAND_MAX;
        x=x*2-1;
        y=(double)rand()/(double)RAND_MAX;
        y=y*2-1;
        distance_squared=x*x+y*y;
        if(distance_squared<=1)
            *num_in_cycle=*num_in_cycle+1;
    }
    MPI_Reduce(num_in_cycle,total_num_in_cycle,1,MPI_LONG_LONG,MPI_SUM,0,comm);
    if(my_rank==0){
        double pi=(double)*total_num_in_cycle/(double)num_point*4;
        printf("the estimate value of pi is %lf\n",pi);
    }
}
