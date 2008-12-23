#include "umtzlib.h"

main()
{
  int i;
  umtzfile* filein;
  umtzfile* fileout;
  float fdata[200];
  filein = umtz_open( "testfile.mtz", "r" );
  fileout = umtz_open( "junk.mtz", "w" );

  printf( "---\n" );
  for (i = 0; i<umtz_num_head(filein); i++)
    printf( "%80.80s\n", umtz_first_head( filein )[i].data );
  printf( "---\n" );
  for (i = 0; i<umtz_num_hist(filein); i++)
    printf( "%80.80s\n", umtz_first_hist( filein )[i].data );
  printf( "---\n" );
  for (i = 0; i<umtz_num_bats(filein); i++)
    printf( "%80.80s\n", umtz_first_bats( filein )[i].cdata );
  printf( "---\n" );

  for (i = 0; i<umtz_num_head(filein); i++)
    umtz_add_head( fileout, umtz_first_head( filein )[i].data );
  for (i = 0; i<umtz_num_hist(filein); i++)
    umtz_add_hist( fileout, umtz_first_hist( filein )[i].data );
  for (i = 0; i<umtz_num_bats(filein); i++)
    umtz_add_bats( fileout, umtz_first_bats( filein )[i].cdata, umtz_first_bats( filein )[i].idata, umtz_first_bats( filein )[i].fdata );

  for (i = 0; i<umtz_num_rows(filein); i++) {
    umtz_get_row(filein, fdata);
    umtz_add_row(fileout, fdata);
    
    if ( i%1000 == 0 ) printf ("%4i %4i %4i %8f %8f\n",(int)fdata[0],(int)fdata[1],(int)fdata[2],fdata[3],fdata[4]);
  }
  umtz_close( filein );
  umtz_close( fileout );
}

