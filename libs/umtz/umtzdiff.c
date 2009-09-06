#include "umtzdiff.h"

main(int argc, char** argv)
{
  mmtzfile file1, file2;
  mmtz_column col; mmtz_dataset set; mmtz_crystal xtl;
  float fdata1[200],fdata2[200],max[200]={0.0},d;
  int i,j,present1[200],present2[200];
  char type[201];

  file1 = mmtz_open( argv[1], "r" );
  file2 = mmtz_open( argv[2], "r" );

  if (mmtz_num_cols(file1) != mmtz_num_cols(file2))
    { for (j = 0; j<mmtz_num_cols(file1); j++) printf("%f ", 999.0); exit(1); }
  if (mmtz_num_rows(file1) != mmtz_num_rows(file2))
    { for (j = 0; j<mmtz_num_cols(file1); j++) printf("%f ", 999.0); exit(1); }

  /* look up column types */
  for (j = 0; j<mmtz_num_cols(file1); j++) {
    mmtz_get_column( file1, j, &col, &set, &xtl );
    type[j]=col.type[0];
  }
  type[j]='#';   /* mark end column */

  /* try and group a phase with another column */
  for (j = 3; j<mmtz_num_cols(file1); j++)
    if ( type[j] == 'P' ) {
      if      ( type[j+1] == 'W' ) { type[j] = 'p'; } /* phi,fom */
      else if ( type[j-1] == 'F' ) { type[j] = 'q'; } /* f, phi */
    }

  /* now read the file data */
  for (i = 0; i<mmtz_num_rows(file1); i++) {
    mmtz_get_row(file1, fdata1, present1);
    mmtz_get_row(file2, fdata2, present2);
    for (j = 0; j<mmtz_num_cols(file1); j++) {
      if ( present1[j] && present2[j] ) {
	switch ( type[j] ) {
	case 'p': /* phi,fom */
	  d = sqrt( pow( fdata1[j+1] * cos( M_PI*fdata1[j]/180.0 ) -
			 fdata2[j+1] * cos( M_PI*fdata2[j]/180.0 ), 2 ) +
		    pow( fdata1[j+1] * sin( M_PI*fdata1[j]/180.0 ) -
			 fdata2[j+1] * sin( M_PI*fdata2[j]/180.0 ), 2 ) );
	  break;
	case 'q': /* f, phi */
	  d = sqrt( pow( fdata1[j-1] * cos( M_PI*fdata1[j]/180.0 ) -
			 fdata2[j-1] * cos( M_PI*fdata2[j]/180.0 ), 2 ) +
		    pow( fdata1[j-1] * sin( M_PI*fdata1[j]/180.0 ) -
			 fdata2[j-1] * sin( M_PI*fdata2[j]/180.0 ), 2 ) );
	  break;
	case 'P': /* any other phase */
	  d = fabs( sin ( M_PI * ( fdata1[j] - fdata2[j] ) / 90.0 ) );
	  break;
	case 'A': /* HL coeff */
	  d = fabs( tanh( fdata1[j] ) - tanh( fdata2[j] ) );
	  break;
	default:
	  d = fabs( fdata1[j] - fdata2[j] );
	}
	if ( d > max[j] ) max[j] = d;
      } else if ( present1[j] ^ present2[j] ) {
	max[j] = 999.0;
      }
    }
  }
  mmtz_close( file1 );
  mmtz_close( file2 );

  for (j = 0; j<mmtz_num_cols(file1); j++) printf( "%f ", max[j] );
}

