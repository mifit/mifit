#include "mmtzlib.h"


mmtzfile mmtz_open( const char* filename, const char* mode )
{ return umtz_open( filename, mode ); }

void mmtz_close( mmtzfile file )
{ umtz_close( file ); }

void mmtz_copy_headers( mmtzfile dest, const mmtzfile src )
{
  umtz_hdr* hpt;
  for ( hpt = umtz_first_head( src ); hpt < umtz_last_head( src ); hpt++ )
    umtz_add_head( dest, hpt->data );
  for ( hpt = umtz_first_hist( src ); hpt < umtz_last_hist( src ); hpt++ )
    umtz_add_hist( dest, hpt->data );
}

int mmtz_num_rows( const mmtzfile file )
{ return umtz_num_rows( file ); }

int mmtz_num_cols( const mmtzfile file )
{ return umtz_num_cols( file ); }

int mmtz_num_datasets( const mmtzfile file )
{ return umtz_num_head_type( file, "DATA" ); }

int mmtz_get_num_symops( const mmtzfile file )
{ return umtz_num_head_type( file, "SYMM" ); }

char* mmtz_get_symop( const mmtzfile file, const int isym, char* symop )
{
  int j = 0;
  umtz_hdr* hpt;
  for ( hpt = umtz_first_head( file ); hpt < umtz_last_head( file ); hpt++ )
    if ( umtz_keymatch( hpt->data, "SYMM" ) ) if ( j++ == isym ) break;
  memcpy( symop, &hpt->data[5], 40 );
  symop[40] = '\0';
  return symop;
}

float* mmtz_get_cell( const mmtzfile file, float* cell )
{
  umtz_get_cell( file, 1, cell );
  return cell;
}

void mmtz_get_setxtl( const mmtzfile file, const int iset, mmtz_dataset* set, mmtz_crystal* xtl )
{
  /* fetch dataset and crystal info */
  int j;
  char ctmp[65];
  float ftmp;
  umtz_hdr* hpt;

  /* initialise info for missing datasets */
  sprintf( set->dname, "unnamed_dataset%i", iset );
  sprintf( xtl->pname, "unnamed_project%i", iset );
  sprintf( xtl->xname, "unnamed_crystal%i", iset );
  set->wavel = 0.0;
  umtz_get_cell( file, iset, xtl->cell );
  for ( hpt = umtz_first_head( file ); hpt < umtz_last_head( file ); hpt++ ) {
    if        ( umtz_keymatch( hpt->data, "CRYS" ) ) {
      sscanf( hpt->data, "%*s %i %64s", &j, ctmp );
      if ( j == iset ) strcpy( xtl->xname, ctmp );
    } else if ( umtz_keymatch( hpt->data, "PROJ" ) ) {
      sscanf( hpt->data, "%*s %i %64s", &j, ctmp );
      if ( j == iset ) strcpy( xtl->pname, ctmp );
    } else if ( umtz_keymatch( hpt->data, "DATA" ) ) {
      sscanf( hpt->data, "%*s %i %64s", &j, ctmp );
      if ( j == iset ) strcpy( set->dname, ctmp );
    } else if ( umtz_keymatch( hpt->data, "DWAV" ) ) {
      sscanf( hpt->data, "%*s %i %f", &j, &ftmp );
      if ( j == iset ) set->wavel = ftmp;
    }
  }
}

int mmtz_get_column( const mmtzfile file, const int icol, mmtz_column* col, mmtz_dataset* set, mmtz_crystal* xtl )
{
  int j = 0, iset = 0;
  umtz_hdr* hpt;

  /* search for the n'th column */
  for ( hpt = umtz_first_head( file ); hpt < umtz_last_head( file ); hpt++ )
    if ( umtz_keymatch( hpt->data, "COLU" ) ) if ( j++ == icol ) break;

  /* grab column info - dataset absent in old files */
  if ( hpt->data[78] == ' ' )
    sscanf( hpt->data, "%*s %30s %2s",
	    col->label, col->type );
  else
    sscanf( hpt->data, "%*s %30s %2s %*f %*f %i",
	    col->label, col->type, &iset );

  /* fetch dataset and crystal info */
   mmtz_get_setxtl( file, iset, set, xtl );

  return icol;
}

int mmtz_add_column( mmtzfile file, const mmtz_column* col, const mmtz_dataset* set, const mmtz_crystal* xtl )
{
  int iset, nset;
  char hdr[200];
  mmtz_dataset fset; mmtz_crystal fxtl;

  /* search for a matching crystal and dataset */
  /* must handle files with DATASET but no CRYSTAL - yeuch */
  nset = mmtz_num_datasets( file );
  for ( iset = 1; iset <= nset; iset++ ) {
    mmtz_get_setxtl( file, iset, &fset, &fxtl );
    if ( ( strcmp( fset.dname, set->dname ) == 0 ) && ( strcmp( fxtl.xname, xtl->xname ) == 0 || strcmp( fxtl.xname, "unnamed_crystal" ) == 0 ) ) break;
  }
  /* if not found, add a new dataset */
  if ( iset > nset ) {
    iset = nset + 1;
    sprintf( hdr, "CRYSTAL  %6i %-64s", iset, xtl->xname );
    umtz_add_head( file, hdr );
    sprintf( hdr, "PROJECT  %6i %-64s", iset, xtl->pname );
    umtz_add_head( file, hdr );
    sprintf( hdr, "DATASET  %6i %-64s", iset, set->dname );
    umtz_add_head( file, hdr );
    sprintf( hdr, "DCELL    %6i %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f", iset, xtl->cell[0], xtl->cell[1], xtl->cell[2], xtl->cell[3], xtl->cell[4], xtl->cell[5] );
    umtz_add_head( file, hdr );
    sprintf( hdr, "DWAVEL   %6i %10.5f", iset, set->wavel );
    umtz_add_head( file, hdr );
  }
  /* add the column */
  sprintf( hdr, "COLUMN %-30s %-2s %16.4f  %16.4f %4i", col->label, col->type, 0.0, 0.0, iset );
  umtz_add_head( file, hdr );

  return ( umtz_num_cols( file ) - 1 );
}

void mmtz_init_headers( mmtzfile file, const char* title, const float cell[6] )
{
  char hdr[80];
  file->headers.size = 0;
  umtz_add_head( file, "VERS MTZ:V1.1 (umtzlib)" );
  sprintf( hdr, "TITLE %-72s", title );
  umtz_add_head( file, hdr );
  sprintf( hdr, "NCOL %8i %12i %8i", 0, 0, 0 );
  umtz_add_head( file, hdr );
  sprintf( hdr, "NDIF %8i", 0 );
  umtz_add_head( file, hdr );
  umtz_add_head( file, "VALM  NAN" );
  sprintf( hdr, "RESO  %11.5f %11.5f ", 0.0, 0.0 );
  umtz_add_head( file, hdr );
  sprintf( hdr, "CELL  %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f", cell[0], cell[1], cell[2], cell[3], cell[4], cell[5] );
  umtz_add_head( file, hdr );
}

void mmtz_add_sort_header( mmtzfile file, const char* sortorder )
{
  char hdr[80];
  int i, s[3];
  for ( i = 0; i < 3; i++ )
    if      ( sortorder[i] == 'H' ) s[i] = 1;
    else if ( sortorder[i] == 'K' ) s[i] = 2;
    else if ( sortorder[i] == 'L' ) s[i] = 3;
  sprintf( hdr, "SORT  %3i %3i %3i %3i %3i", s[0], s[1], s[2], 0, 0 );
  umtz_add_head( file, hdr );
}

void mmtz_add_syminf_header( mmtzfile file, const int nsym, const int nsymp, const char cellcode, const int ccp4_symbol, const char* hm_symbol, const char* pg_symbol )
{
  char hdr[80];
  sprintf( hdr, "SYMINF %3i %2i %c %5i %-10s %s", nsym, nsymp, cellcode, ccp4_symbol, hm_symbol, pg_symbol );
  umtz_add_head( file, hdr );
}

void mmtz_add_symop_header( mmtzfile file, const char* symop )
{
  char hdr[80];
  sprintf( hdr, "SYMM %s", symop );
  umtz_add_head( file, hdr );
}

void mmtz_get_row( const mmtzfile file, float* fdata, int* flags )
{
  int i;
  umtz_get_row( file, fdata );
  for ( i = 0; i < umtz_num_cols( file ); i++ )
    flags[i] = !umtz_ismnf( file, fdata[i] );
}

void mmtz_add_row( mmtzfile file, float* fdata, const int* flags )
{
  int i;
  for ( i = 0; i < umtz_num_cols( file ); i++ )
    if ( !flags[i] ) fdata[i] = umtz_mnf( file );
  umtz_add_row( file, fdata );
}

void mmtz_seek_row( const mmtzfile file, const int n )
{
  umtz_seek_row( file, n );
}
