#include <mpi.h>                                        //line 1
#include <stdio.h>                                      //line 2
int main( int argc, char **argv ) {                     //line 3
    
    int global_p = 0;
    
    int rank, size;                                     //line 4
    MPI_Init( &argc, &argv );                           //line 5
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );             //line 6
    MPI_Comm_size( MPI_COMM_WORLD, &size );             //line 7
    
    int Pavg = 0;
    int Tavg = 0;
    int i;
    int novaVar;
    
    if(rank==0){
        printf("QUANTAS VEZES É QUE VOU APARECER NO ECRÃ??\n");
        novaVar = 50;
        //printf("Novavar: %d\n", novaVar);
    }
    else{
        printf(",,,");
    }

    if(rank==0){
        printf( "Hello from root process %d/%d\n", rank, size ); //line 8
        Pavg += 2;
        i = 10;
        printf("Pavg: %d\n", Pavg);
        printf("Novavar: %d\n", novaVar);
        //global_p += 1;
    }
    else{
        printf( "Hello from not root process %d/%d\n", rank, size ); //line 8
        Pavg += 1;
        printf("Pavg: %d\n", Pavg);
    }
    
    printf("i: %d\n", i);

    MPI_Barrier(MPI_COMM_WORLD);

    printf("before global_p: %d, on process %d\n", global_p, rank);

    MPI_Reduce(&Pavg, &global_p, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    MPI_Finalize();                                    //line 9
    printf("after global_p: %d, on process %d\n", global_p, rank);
    return 0;                                           //line 10
}