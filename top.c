#include <mpi.h>
#include <stdio.h>
MPI_Comm  old_comm, new_comm;
int ndims, reorder, ierr;
int dim_size[2], periods[2];

int main(int argc, char *argv[])
{
  int rank;
  MPI_Comm new_comm,linecomm;	
  MPI_Init(&argc,&argv);	
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  old_comm = MPI_COMM_WORLD;
  ndims = 2;
  dim_size[0] = 3;
  dim_size[1] = 2;
  periods[0] = 1;
  periods[1] = 0;
  reorder = 1;
      
  ierr =  MPI_Cart_create(old_comm,ndims,dim_size,
             periods,reorder,&new_comm);

  int coords[2];
  int rankx, ranky;

  MPI_Cart_get(new_comm,2,dim_size,periods,coords);
  ranky=coords[0]; rankx=coords[1];

  printf("global %d cart x %d y %d \n",rank,rankx,ranky);

//  int free_coords[] = {0,1};
//  MPI_Cart_sub(new_comm,free_coords,linecomm);



   MPI_Finalize();
   
  return 0;
}
