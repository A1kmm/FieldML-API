/******************************************************************************
 **  PARALLEL HDF5 TEST SUITE
 **
 **  C++ suite that tests the parallel HDF5 APIs added to FieldML library.
 **
 **  Mark Cheeseman, NIWA
 **  January 14, 2013
 ******************************************************************************/

#include <mpi.h>
#include <math.h>
#include "hdf5.h"
#include "FieldmlIoApi.h"
#include "fieldml_api.h" 

const int global_arraysize = 68921;
const int global_lowres_arraysize = 1331;
const int array_length = 41;
const int lowres_array_length = 11;
const int num_dim = 3;
const int num_time_slices = 100;

const double maximum_time = 40.0;
const double propagation_speed = maximum_time / (double)( num_time_slices - 1 );
const double radius = ((double)(array_length-1))/2.0;
const double dt = maximum_time / (double)(num_time_slices-1);

int count[2], offset[2], low_count[2], low_offset[2];
double **highres_coordinates, **lowres_coordinates, **p_highres;


/*=============================================================================
  DECOMPOSE_DOMAIN
 
  C++ subroutine that determines the array size on each MPI task 
  =============================================================================*/
 
void decompose_domain( int myID, int num_mpi_tasks ) {

     int remainder;

     count[1] = global_arraysize / num_mpi_tasks;
     offset[0] = 0;
     offset[1] = count[1]*myID;
     remainder = global_arraysize % num_mpi_tasks;
     if ( remainder>0 ) {
        if ( myID<remainder ) { count[1]++; }
        else { offset[1] += remainder-1; }
     }

     low_count[1] = global_lowres_arraysize / num_mpi_tasks;
     low_offset[0] = 0;
     low_offset[1] = low_count[1]*myID;
     remainder = global_arraysize % num_mpi_tasks;
     if ( remainder>0 ) {
        if ( myID<remainder ) { low_count[1]++; }
        else { low_offset[1] += remainder-1; }
     }
     return;

}

#include "coordinate.cpp"

/*=============================================================================
  COMPUTE_POTENTIAL
 
  C++ subroutine that computes the value of the potential field at a given 
  point (x,y,z) and time (t). 
  =============================================================================*/

double compute_potential( double x, double y, double z, double t ) {

       double x_dist, y_dist, z_dist, distance, potential;

       x_dist = x - radius;
       y_dist = y - radius;
       z_dist = z - radius;
       distance = sqrt( x_dist*x_dist + y_dist*y_dist + z_dist*z_dist );
       potential = cos(  (1.570796325 * ( distance - (propagation_speed * t)) / radius ) );
       if ( potential<0.0 ) { potential = 0.0; }
       return potential;

}

void construct_highres_potential() {

     int i,j;
     double time = 0.0;

     p_highres = new double*[num_time_slices];
     for ( i=0; i<num_time_slices; i++ ) 
         p_highres[i] = new double[count[1]]; 
     
     for ( i=0; i<num_time_slices; i++ ) {
         for ( j=0; j<count[1]; j++ ) {
             p_highres[i][j] = compute_potential( highres_coordinates[0][j],
                                                  highres_coordinates[1][j],
                                                  highres_coordinates[2][j],
                                                  time );
         }
         time += dt;
     }

}

/*=============================================================================
  WRITE_POTENTIAL_DATA
 
  C++ subroutine that writes the potential data array to hard disk.
  =============================================================================*/

bool write_potential_data() {

     bool testOK;
     int tmp;
     FmlSessionHandle session;
     FmlObjectHandle cType, resource, sourceD, writer;
     FmlIoErrorNumber status;

 /*
  * Create a HDF5 file for the coordinate data and ready it for parallel IO
  * access.
  *----------------------------------------------------------------------------*/
     session = Fieldml_Create( "test", "test" );
     cType = Fieldml_CreateContinuousType( session, "test.scalar_real" );

     resource = Fieldml_CreateHrefDataResource( session, "potential.resource",
                                               "PHDF5", "potential.h5" );
     sourceD = Fieldml_CreateArrayDataSource( session, "potential.source_double",
                                              resource, "highres_potential", 2 );

 /*
  * Reopen HDF5 file and write the high-resolution data into it
  *----------------------------------------------------------------------------*/
     tmp = count[1];
     count[0] = num_time_slices;
     count[1] = global_arraysize;
     writer = Fieldml_OpenArrayWriter( session, sourceD, cType, 0, count, 2 );

     count[1] = tmp;
     status = Fieldml_WriteDoubleSlab( writer, offset, count, &(p_highres[0][0]) );
     if ( status==FML_IOERR_NO_ERROR ) { testOK = true; }
     else { testOK = false; }

     Fieldml_CloseWriter( writer );
     Fieldml_Destroy( session );
     return testOK;
}

/**                                                                          **
 **>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MAIN PROGRAM >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>**
 **                                                                          **/

int main( int argc, char *argv[] ) {

    int num_mpi_tasks, myID;
    bool status;
    
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &num_mpi_tasks );
    MPI_Comm_rank( MPI_COMM_WORLD, &myID );
    
    if ( myID==0 ) {
       printf( "\n------------------------------------------------------------\n" );
       printf( "      PARALLEL HDF5 TEST SUITE v1.0\n" );
       printf( "------------------------------------------------------------\n\n" );
    }

    decompose_domain( myID, num_mpi_tasks );
   
/***---------------------------------------------------------------------------
 *** HIGH RESOLUTION TESTS
 ***---------------------------------------------------------------------------*/
    
    construct_coordinate_data();
    construct_highres_potential();

    status = write_coordinate_data();
    if ( myID==0 ) {
       if ( status ) { printf( "   High Resolution Tests:  coordinates... SUCCESS\n" ); }
       else { printf( "   High Resolution Tests:  coordinates... FAIL\n" ); }
    }

    status = write_potential_data();
    if ( myID==0 ) {
       if ( status ) { printf( "                           potential..... SUCCESS\n\n" ); }
       else { printf( "                           potential..... FAIL\n\n" ); }
    }



    MPI_Finalize();
    return 0;
}
