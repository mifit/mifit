#include "mmtzlib.h"

main()
{
  int i;
  float fdata[200];
  int idata[200];
  char cdata[80];
  mmtzfile filein;
  mmtzfile fileout;
  mmtz_column col;
  mmtz_crystal xtl;
  mmtz_dataset set;
  filein  = mmtz_open( "testfile.mtz", "r" );
  fileout = mmtz_open( "junk.mtz", "w" );

  for ( i = 0; i < mmtz_num_cols( filein ); i++ ) {
    mmtz_get_column( filein, i, &col, &set, &xtl );
    printf("%2i %16.16s %16.16s %16.16s %8.3f\n", i, col.label, set.dname, xtl.xname, set.wavel);
  }
  for ( i = 0; i < mmtz_get_num_symops( filein ); i++ ) {
    mmtz_get_symop( filein, i, cdata );
    printf("%2i %s\n",i,cdata);
  }

  mmtz_copy_headers( fileout, filein );

  strcpy( col.label, "newcol-oldset" );
  mmtz_add_column( fileout, &col, &set, &xtl );
  strcpy( col.label, "newcol-newset" );
  strcpy( set.dname, "newset" ); set.wavel = 9.9;
  mmtz_add_column( fileout, &col, &set, &xtl );

  for (i = 0; i<mmtz_num_rows(filein); i++) {
    mmtz_get_row(filein, fdata, idata);
    mmtz_add_row(fileout, fdata, idata);
  }
  mmtz_close( filein );
  mmtz_close( fileout );
}
