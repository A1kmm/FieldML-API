/**************************************************************************************
 *** SIMPLE_HDF5_TEST
 ***
 *** C program that tests the serial HDF5 library on a standard FieldML file.
 ***
 ***   Mark Cheeseman, NIWA
 ***   January 10, 2012
 **************************************************************************************/

#include <stdlib.h>
#include <unistd.h>
#include "hdf5.h"

int main( int argc, char **argv ) {

    int n;
    hid_t fid, dset;
    herr_t status;
    short iarray[20];
    double darray[12];

    status = 0;
    if ( argc>1 ) {
       if ( access(argv[1],R_OK) != -1 ) {
          printf( "File to be tested: %s\n\n", argv[1] );
          fid = H5Fopen( argv[1], H5F_ACC_RDONLY, H5P_DEFAULT );
          dset = H5Dopen( fid, "/I16BE", H5P_DEFAULT );
          status = H5Dread( dset, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iarray ); 
          if ( status<0 ) { 
             printf( "ERROR - read operation failed\n" ); 
             exit(1);
          }
          status = H5Dclose( dset );
          dset = H5Dopen( fid, "/DOUBLE", H5P_DEFAULT );
          status = H5Dread( dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, darray ); 
          status = H5Dclose( dset );
          status = H5Fclose( fid );

          printf( "values in I16BE were:\n" );
          for ( n=0; n<20; n++ ) { printf( "%d\n", iarray[n] ); }
          printf( "values in DOUBLE were:\n" );
          for ( n=0; n<12; n++ ) { printf( "%7.3f\n", darray[n] ); }
       }
    } else {
       printf( "ERROR - an input filename must be specified!\n\n" );
       exit(1);
    }      
    return 0;
}
