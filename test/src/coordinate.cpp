/*=============================================================================
  CONSTRUCT_COORDINATE_DATA
 
  C++ subroutine that forms the coordinate data array 
  =============================================================================*/

void construct_coordinate_data() {
     
     int i,j,k,global_id,local_id;
     double tmp, zValue, yValue, xValue;

 /*
  * Create the high and low-resolution coordinate arrays
  *----------------------------------------------------------------------------*/
     highres_coordinates = new double*[num_dim];
     lowres_coordinates = new double*[num_dim];

     for ( i=0; i<num_dim; i++ ) { 
         highres_coordinates[i] = new double[count[1]]; 
         lowres_coordinates[i] = new double[count[1]]; 
     }

 /*
  * Construct the high-resolution data values
  *----------------------------------------------------------------------------*/
     tmp = (double)lowres_array_length / ((double)lowres_array_length - 1.0);
     global_id = 0;
     for ( k=0; k<lowres_array_length; k++ ) {
         zValue = tmp * (double)k;
         for ( j=0; j<lowres_array_length; j++ ) {
             yValue = tmp * (double)j;
             for ( i=0; i<lowres_array_length; i++ ) {
                 xValue = tmp * (double)i;
                 if ( (global_id>=low_offset[1])&&(global_id<low_offset[1]+low_count[1]) ) {
                    local_id = global_id - low_offset[1]; 
                    lowres_coordinates[0][local_id] = xValue; 
                    lowres_coordinates[1][local_id] = yValue; 
                    lowres_coordinates[2][local_id] = zValue; 
                 }
                 global_id++;
             }
         }
     }

 /*
  * Construct the low-resolution data values
  *----------------------------------------------------------------------------*/
     tmp = (double)array_length / ((double)array_length - 1.0);
     global_id = 0;
     for ( k=0; k<array_length; k++ ) {
         zValue = tmp * (double)k;
         for ( j=0; j<array_length; j++ ) {
             yValue = tmp * (double)j;
             for ( i=0; i<array_length; i++ ) {
                 xValue = tmp * (double)i;
                 if ( (global_id>=offset[1])&&(global_id<offset[1]+count[1]) ) {
                    local_id = global_id - offset[1]; 
                    highres_coordinates[0][local_id] = xValue; 
                    highres_coordinates[1][local_id] = yValue; 
                    highres_coordinates[2][local_id] = zValue; 
                 }
                 global_id++;
             }
         }
     }

     return;

}

/*=============================================================================
  WRITE_COORDINATE_DATA
 
  C++ subroutine that writes the coordinate data array to hard disk.
  =============================================================================*/

bool write_coordinate_data() {

     bool testOK;
     int tmp;
     FmlSessionHandle session;
     FmlObjectHandle cType, resource, sourceD, sourceD2, writer;
     FmlIoErrorNumber status;

 /*
  * Create a HDF5 file for the coordinate data and ready it for parallel IO
  * access.
  *----------------------------------------------------------------------------*/
     session = Fieldml_Create( "test", "test" );
     cType = Fieldml_CreateContinuousType( session, "test.scalar_real" );

     resource = Fieldml_CreateHrefDataResource( session, "coordinates.resource",
                                               "PHDF5", "coordinates.h5" );
     sourceD = Fieldml_CreateArrayDataSource( session, "coordinates.source_double",
                                              resource, "highres_coordinates", 2 );

 /*
  * Reopen HDF5 file and write the high-resolution data into it
  *----------------------------------------------------------------------------*/
     tmp = count[1];
     count[0] = num_dim;
     count[1] = global_arraysize;
     writer = Fieldml_OpenArrayWriter( session, sourceD, cType, 0, count, 2 );

     count[1] = tmp;
     status = Fieldml_WriteDoubleSlab( writer, offset, count, &(highres_coordinates[0][0]) );
     if ( status==FML_IOERR_NO_ERROR ) { testOK = true; }
     else { testOK = false; }

     Fieldml_CloseWriter( writer );
     Fieldml_Destroy( session );
     return testOK;
}
