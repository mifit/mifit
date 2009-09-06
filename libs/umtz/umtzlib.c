#include "umtzlib.h"

int umtz_num_bats( const umtzfile* file );

umtzfile* umtz_open( const char* filename, const char* cmode )
{
  int i, j, nhst, nint, nflt, arch, hdr_off;
  umtzfile* file;
  char hdr[240]; int idata[MTZBATINT]; float fdata[MTZBATFLT];

  file = (umtzfile*)malloc( sizeof(umtzfile) );
  file->ncol = file->nrow = 0;
  file->mnf = ccp4_nan().f;
  file->headers.size = file->history.size = file->batches.size = 0;
  file->headers.capacity = file->history.capacity = file->batches.capacity = 4;
  file->headers.data = malloc( file->headers.capacity * sizeof(umtz_hdr) );
  file->history.data = malloc( file->history.capacity * sizeof(umtz_hdr) );
  file->batches.data = malloc( file->batches.capacity * sizeof(umtz_bat) );
  strncpy( file->mode, cmode, 4 );
  strncpy( file->filename, filename, 200 );

  if ( file->mode[0] == 'r' ) {
    /*
      OPEN AN MTZ FOR READING
      Read and store headers
    */
    file->iunit = ccp4_qopen( filename, 3 );
    /* reach arch */
    arch = ccp4_qrarch( file->iunit, 2 );
    ccp4_qseek( file->iunit, 0, 0, 1 );
    /* read stamp */
    ccp4_qmode( file->iunit, 0 );
    ccp4_qreadc( file->iunit, hdr, 4 );
    if ( !umtz_keymatch( hdr, "MTZ" ) ) CCP4_fatal("Not an mtz");
    /* read header offset */
    ccp4_qmode( file->iunit, 6 );
    ccp4_qread( file->iunit, (uint8 *) &hdr_off, 1 );
    ccp4_qseek( file->iunit, 0, hdr_off-1, 1 );
    ccp4_qmode( file->iunit, 0 );
    /* read headers */
    i = ccp4_qreadc( file->iunit, hdr, MTZRECLEN );
    while ( i >= 0 && !umtz_keymatch( hdr, "END" ) ) {
      umtz_add_head( file, hdr );
      i = ccp4_qreadc( file->iunit, hdr, MTZRECLEN );
    }
    /* read optional blocks */
    while ( i >= 0 && !umtz_keymatch( hdr, "MTZE") ) {
      if ( umtz_keymatch( hdr, "MTZHIST" ) ) {
	/* read optional history block */
	sscanf( hdr, "MTZHIST %i", &nhst );
	for ( j = 0; j < nhst; j++ ) {
	  ccp4_qreadc( file->iunit, hdr, MTZRECLEN );
	  umtz_add_hist( file, hdr ); /* read next line */
	}
	i = ccp4_qreadc( file->iunit, hdr, MTZRECLEN );
      } else if ( umtz_keymatch( hdr, "MTZBATS") ) {
	/* read optional batch block */
	i = ccp4_qreadc( file->iunit, hdr, MTZRECLEN );
	while ( i >= 0 && umtz_keymatch( hdr, "BH") ) {
	  sscanf( hdr, "BH %*i %*i %i %i", &nint, &nflt );
	  i = ccp4_qreadc( file->iunit, &hdr[MTZRECLEN], MTZRECLEN );
	  ccp4_qmode( file->iunit, 6 );
	  i = ccp4_qread( file->iunit, (uint8 *)idata, nint);
	  ccp4_qmode( file->iunit, 2 );
	  i = ccp4_qread( file->iunit, (uint8 *)fdata, nflt);
	  ccp4_qmode( file->iunit, 0 );
	  i = ccp4_qreadc( file->iunit, &hdr[2*MTZRECLEN], MTZRECLEN );
	  umtz_add_bats( file, hdr, idata, fdata );
	  i = ccp4_qreadc( file->iunit, hdr, MTZRECLEN ); /* read next line */
	}
      } else {
	i = ccp4_qreadc( file->iunit, hdr, MTZRECLEN ); /* read next line */
      }
    }
    /* calc number of reflections */
    file->nrow = (hdr_off - MTZDATAOFF - 1) / umtz_num_cols( file );

  } else if ( file->mode[0] == 'w' ) {
    /*
      OPEN AN MTZ FOR WRITING
      Write null header block
    */
    file->iunit = ccp4_qopen( filename, 4 );
    /* note we will write the header info later */
  } else CCP4_fatal("umtzlib: unknown file mode: must be 'r' or 'w'");

  /* set up for reading/writing reflections */
  ccp4_qmode( file->iunit, 2 );
  ccp4_qseek( file->iunit, 0, MTZDATAOFF, 1 );
  return file;
}

void umtz_close( umtzfile* file )
{
  int nint, nflt, hdr_off;
  char hdr[81];
  umtz_hdr* hpt; umtz_bat* bpt;

  if ( file->mode[0] == 'r' ) {
    /*
      CLOSE AN MTZ FOR READING
    */
    ccp4_qclose( file->iunit );

  } else {
    /*
      CLOSE AN MTZ FOR WRITING
      Fill in both initial and terminal headers
    */

    /* write stamp */
    ccp4_qseek( file->iunit, 0, 0, 1 );
    ccp4_qmode( file->iunit, 0 );
    ccp4_qwritc( file->iunit, "MTZ ", 4 );
    /* header offset */
    hdr_off = file->nrow * file->ncol + MTZDATAOFF + 1;
    ccp4_qmode( file->iunit, 2 );
    ccp4_qwrite( file->iunit, (uint8 *) &hdr_off, 1 );
    /* write arch */
    ccp4_qwarch( file->iunit, 2 );

    /* rewrite column min/max information */
    umtz_rewrite_headers_ranges( file );
    /* rewrite redundant legacy headers */
    umtz_rewrite_headers_legacy( file );

    /* write terminal headers */
    ccp4_qseek( file->iunit, 0, hdr_off-1, 1 );
    ccp4_qmode( file->iunit, 0 );
    for ( hpt = umtz_first_head( file ); hpt < umtz_last_head( file ); hpt++ )
      ccp4_qwritc( file->iunit, hpt->data, MTZRECLEN );
    ccp4_qwritc( file->iunit, "END                                                                             ", MTZRECLEN );
    sprintf( hdr, "MTZHIST %3i                                                                     ", umtz_num_hist( file ) );
    ccp4_qwritc( file->iunit, hdr , MTZRECLEN );
    for ( hpt = umtz_first_hist( file ); hpt < umtz_last_hist( file ); hpt++ )
      ccp4_qwritc( file->iunit, hpt->data, MTZRECLEN );
    if ( umtz_num_bats( file ) > 0 ) ccp4_qwritc( file->iunit, "MTZBATS                                                                         ", MTZRECLEN );
    for ( bpt = umtz_first_bats( file ); bpt < umtz_last_bats( file ); bpt++ ) {
      sscanf( bpt->cdata, "BH %*i %*i %i %i", &nint, &nflt );
      ccp4_qwritc( file->iunit, &bpt->cdata[0], MTZRECLEN );
      ccp4_qwritc( file->iunit, &bpt->cdata[MTZRECLEN], MTZRECLEN );
      ccp4_qmode( file->iunit, 6 );
      ccp4_qwrite( file->iunit, (uint8 *)bpt->idata, nint );
      ccp4_qmode( file->iunit, 2 );
      ccp4_qwrite( file->iunit, (uint8 *)bpt->fdata, nflt );
      ccp4_qmode( file->iunit, 0 );
      ccp4_qwritc( file->iunit, &bpt->cdata[2*MTZRECLEN], MTZRECLEN );
    }
    ccp4_qwritc( file->iunit, "MTZENDOFHEADERS                                                                 ", MTZRECLEN );
    ccp4_qclose( file->iunit );
  }

  /* tidy up */
  free( file->headers.data );
  free( file->history.data );
  free( file->batches.data );
  free( file );
}

int umtz_num_cols( const umtzfile* file ) { return file->ncol; }

int umtz_num_rows( const umtzfile* file ) { return file->nrow; }

int umtz_num_head( const umtzfile* file ) { return file->headers.size; }

int umtz_num_hist( const umtzfile* file ) { return file->history.size; }

int umtz_num_bats( const umtzfile* file ) { return file->batches.size; }

float umtz_mnf( const umtzfile* file ) { return file->mnf; }

int umtz_ismnf( const umtzfile* file, float f )
{
  union float_uint_uchar u;
  u.f = file->mnf;
  if ( ccp4_isnan( &u ) ) { u.f = f; return ccp4_isnan( &u ); }
  return ( f == file->mnf );
}

umtz_hdr* umtz_first_head( const umtzfile* file )
{ return (umtz_hdr*)file->headers.data; }

umtz_hdr* umtz_first_hist( const umtzfile* file )
{ return (umtz_hdr*)file->history.data; }

umtz_bat* umtz_first_bats( const umtzfile* file )
{ return (umtz_bat*)file->batches.data; }

umtz_hdr* umtz_last_head( const umtzfile* file )
{ return (umtz_hdr*)file->headers.data + file->headers.size; }

umtz_hdr* umtz_last_hist( const umtzfile* file )
{ return (umtz_hdr*)file->history.data + file->history.size; }

umtz_bat* umtz_last_bats( const umtzfile* file )
{ return (umtz_bat*)file->batches.data + file->batches.size; }

int umtz_keymatch( const char* hdr, const char* key )
{ return ( strncmp( hdr, key, strlen(key) ) == 0 ); }

void umtz_add_head( umtzfile* file, const char* hdr )
{
  /* cant add to a header once file has reflections (can't change ncol) */
  if ( file->nrow != 0 ) CCP4_fatal("Can't add headers after reflections");
  /* if we add a column keyword, number of columns increases */
  if ( umtz_keymatch( hdr, "COLU" ) ) file->ncol++;
  /* if we add a VALM keyword, update mnf. (non-string leaves as NaN) */
  if ( umtz_keymatch( hdr, "VALM" ) ) sscanf( hdr, "VALM %f", &file->mnf );
  /* add the column */
  umtz_make_rec( &file->headers, MTZRECLEN );
  umtz_copy_pad( umtz_first_head( file )[umtz_num_head(file)-1].data, hdr, MTZRECLEN );
}

void umtz_add_hist( umtzfile* file, const char* hdr )
{
  umtz_make_rec( &file->history, MTZRECLEN );
  umtz_copy_pad( umtz_first_hist( file )[umtz_num_hist(file)-1].data, hdr, MTZRECLEN );
}

void umtz_add_bats( umtzfile* file, const char* cdata, const int* idata, const float* fdata )
{
  int nint, nflt;
  sscanf( cdata, "BH %*i %*i %i %i", &nint, &nflt );
  umtz_make_rec( &file->batches, MTZBATLEN );
  memcpy( umtz_first_bats( file )[umtz_num_bats(file)-1].cdata, cdata,MTZBATCHR);
  memcpy( umtz_first_bats( file )[umtz_num_bats(file)-1].idata, idata, 4*nint );
  memcpy( umtz_first_bats( file )[umtz_num_bats(file)-1].fdata, fdata, 4*nflt );
}

void umtz_seek_row( const umtzfile* file, const int n )
{ ccp4_qseek( file->iunit, 0, MTZDATAOFF + n * file->ncol, 1 ); }

void umtz_get_row( const umtzfile* file, float* fdata )
{ ccp4_qread( file->iunit, (uint8 *) fdata, file->ncol ); }

void umtz_add_row( umtzfile* file, const float* fdata )
{
  if ( file->mode[0] != 'w' ) CCP4_fatal("umtzlib: cannot write to read file");
  ccp4_qseek( file->iunit, 0, MTZDATAOFF + file->nrow * file->ncol, 1 );
  ccp4_qwrite( file->iunit, (uint8 *) fdata, file->ncol );
  file->nrow++;
}

void umtz_get_cell( const umtzfile* file, const int ixtl, float* cell )
{
  /* get the base cell from the file */
  /* FIXME: currently gets the cell of the first dataset, or if that
     does not exist, the CELL line. Is that right? */
  int j;
  umtz_hdr* hpt;

  for ( hpt = umtz_first_head( file ); hpt < umtz_last_head( file ); hpt++ ) {
    if        ( umtz_keymatch( hpt->data, "CELL" ) ) {
      sscanf( hpt->data, "CELL %f %f %f %f %f %f", &cell[0], &cell[1], &cell[2], &cell[3], &cell[4], &cell[5] );
      if ( ixtl == 0 ) break;
    } else if ( umtz_keymatch( hpt->data, "DCEL" ) ) {
      sscanf( hpt->data, "DCELL %i %f %f %f %f %f %f", &j, &cell[0], &cell[1], &cell[2], &cell[3], &cell[4], &cell[5] );
      if ( j == ixtl ) break;
    }
  }
}

void umtz_cell_metric( const float* cell, float* metric )
{
  double a, b, c, alph, beta, gamm, as, bs, cs, vol;
  /* calc metric tensor */
  a = cell[0]; alph = cell[3] * M_PI / 180.0;
  b = cell[1]; beta = cell[4] * M_PI / 180.0;
  c = cell[2]; gamm = cell[5] * M_PI / 180.0;
  vol = a*b*c*sqrt( 1.0 - cos(alph)*cos(alph) - cos(beta)*cos(beta) - cos(gamm)*cos(gamm) + 2.0*cos(alph)*cos(beta)*cos(gamm) );
  as = b*c*sin(alph)/vol;
  bs = c*a*sin(beta)/vol;
  cs = a*b*sin(gamm)/vol;
  metric[0] = (float)(as*as);
  metric[1] = (float)(bs*bs);
  metric[2] = (float)(cs*cs);
  metric[3] = (float)(2.0*as*bs*(cos(beta)*cos(alph)-cos(gamm)) / (sin(alph)*sin(beta)));
  metric[4] = (float)(2.0*cs*as*(cos(alph)*cos(gamm)-cos(beta)) / (sin(gamm)*sin(alph)));
  metric[5] = (float)(2.0*bs*cs*(cos(gamm)*cos(beta)-cos(alph)) / (sin(beta)*sin(gamm)));
}

int umtz_num_head_type( const umtzfile* file, char* head )
{
  umtz_hdr* hpt;
  int nhead = 0;
  for ( hpt = umtz_first_head( file ); hpt < umtz_last_head( file ); hpt++ )
    if ( umtz_keymatch( hpt->data, head ) ) nhead++;
  return nhead;
}

void umtz_make_rec( umtz_list* l, const int rec_len )
{
  if ( l->size == l->capacity ) {
    l->capacity = 2 * l->capacity;
    l->data = realloc( l->data, l->capacity * rec_len );
  }
  l->size++;
}

void umtz_copy_pad( char* d, const char* s, const int rec_len  )
{
  int i;
  for ( i = 0; i < rec_len; i++ ) { d[i] = s[i]; if (d[i]=='\0') break; }
  for (      ; i < rec_len; i++ ) { d[i] = ' '; }
}

void umtz_rewrite_headers_ranges( umtzfile* file )
{
  /* Do any programs use the column min/max information? If not, the
     information should be deprected and this code removed */
  int i, j;
  double s, smax=0.0, smin=1.0e8;
  float* min; float* max; float* fdata; float cell[6], metric[6];
  char txt[40];
  umtz_hdr* hpt;

  min = (float*)malloc( file->ncol * sizeof(float) );
  max = (float*)malloc( file->ncol * sizeof(float) );
  fdata = (float*)malloc( file->ncol * sizeof(float) );
  for ( i = 0; i < file->ncol; i++ ) { min[i] = 999.0; max[i] = -999.0; }

  /* get cell */
  umtz_get_cell( file, 0, cell );
  umtz_cell_metric( cell, metric );

  /* read the reflection data */
  ccp4_qseek( file->iunit, 0, MTZDATAOFF, 1 );
  ccp4_qmode( file->iunit, 2 );
  for ( j = 0; j < file->nrow; j++ ) {
    umtz_get_row( file, fdata );
    for ( i = 0; i < file->ncol; i++ ) {
      if ( !umtz_ismnf( file, fdata[i] ) ) {
	if ( fdata[i] < min[i] ) min[i] = fdata[i]; /* min and max values */
	if ( fdata[i] > max[i] ) max[i] = fdata[i];
      }
    }
    s = metric[0]*fdata[0]*fdata[0] + metric[1]*fdata[1]*fdata[1]
      + metric[2]*fdata[2]*fdata[2] + metric[3]*fdata[0]*fdata[1]
      + metric[4]*fdata[0]*fdata[2] + metric[5]*fdata[1]*fdata[2];
    if ( s < smin ) smin = s; /* min and max resolution */
    if ( s > smax ) smax = s;
  }

  /* now update the headers */
  j = 0; /* j is column number */
  for ( hpt = umtz_first_head( file ); hpt < umtz_last_head( file ); hpt++ ) {
    if ( umtz_keymatch( hpt->data, "RESO" ) ) {
      sprintf( txt, "RESO  %11.5f %11.5f ", smin, smax );
      memcpy( hpt->data, txt, 29 );
    }
    if ( umtz_keymatch( hpt->data, "COLU" ) ) {
      sprintf( txt, " %16.4f  %16.4f", min[j], max[j] );
      memcpy( hpt->data + 40, txt, 35 );
      j++;
    }
  }

  free( min );
  free( max );
  free( fdata );
}

void umtz_rewrite_headers_legacy( umtzfile* file )
{
  /* The mtz headers contain redundant legacy sizing info, which was
     needed before dynamic memory became available. umtz does not read
     these headers, and only writes them for compatibility with other
     libraries. These headers should be deprecated. */
  /*
    Redundant headers include:
     NCOL: number of columns and batches are determined from the
           number of COLUMN and BATCH headers. nref = (hdr_off -
           MTZDATAOFF - 1)/ncol
     NDIF: number of datasets determined from dataset headers
     MTZH: number of history entries determined by counting them
  */
  int nbat, nset;
  char txt[80];
  umtz_hdr* hpt;

  /* count batches, datasets, history */
  nset = umtz_num_head_type( file, "DATA" );
  nbat = umtz_num_bats( file );

  /* rewrite headers */
  for ( hpt = umtz_first_head( file ); hpt < umtz_last_head( file ); hpt++ ) {
    if ( umtz_keymatch( hpt->data, "NCOL" ) ) {
      sprintf( txt, "NCOL %8i %12i %8i ", file->ncol, file->nrow, nbat );
      memcpy( hpt->data, txt, 35 );
    }
    if ( umtz_keymatch( hpt->data, "NDIF" ) ) {
      sprintf( txt, "NDIF %8i   ", nset );
      memcpy( hpt->data, txt, 15 );
    }
  }
}
