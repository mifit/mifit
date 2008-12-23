#if 0
/*! \file mmtzlib.h \brief The mini-mtz user-level library.

This is a simple lightweight library which can be used by application
developers to access MTZ files from C, and other languages which can
access C APIs.

Accessing a file is achieved by use of an mmtzfile structure. This
structure is allocated by calling the mmtz_open() function with a
filename and either "r" or "w" for read or write access. Closing the
file frees the structure.

When reading a file, this structure is the queried to obtain
information from the MTZ headers. The actual data in the file may be
accessed by using the mmtz_get_row() function to access the data for
each reflection in turn.

When writing a file, the file is opened, and the headers must be
constructed using the 'init' and 'add' functions. Only once
construction of the headers is complete may any actual data be written
to the file. Once all the data is written, the file must be closed to
force the headers to be written.

To append columns to a file, open on file for reading and a second for
writing. Use the mmtz_copy_headers() function to copy the headers from
the first to the second. Then add any additional column headers using
mmtz_add_column(). Finally read the data a reflection at a time, fill
in the additional columns and write the data. The resulting code looks
like this:
\code
  mmtzfile filein, fileout;
  mmtz_column col;
  mmtz_crystal xtl;
  mmtz_dataset set;
  int fcol, pcol;

  filein  = mmtz_open( "old.mtz", "r" );
  fileout = mmtz_open( "new.mtz", "w" );
  mmtz_copy_headers( fileout, filein );
  \/\* get dataset and crystal for column 0 (i.e. H) \*\/
  mmtz_get_column( filein, 1, &col, &set, &xtl );
  strcpy( col.label, "FCAL" );    \/\* add an FCAL column \*\/
  strcpy( col.type, "F" );
  fcol = mmtz_add_column( fileout, &col, &set, &xtl );
  strcpy( col.label, "PHICAL" );  \/\* add an PHICAL column \*\/
  strcpy( col.type, "P" );
  pcol = mmtz_add_column( fileout, &col, &set, &xtl );
  for (i = 0; i<mmtz_num_rows(filein); i++) {
    mmtz_get_row(filein, fdata, idata);
    fdata[ fcol ]   = 500.0;
    fdata[ pcol ] = 60.0;
    mmtz_add_row(fileout, fdata, idata);
  }
  mmtz_close( filein );
  mmtz_close( fileout );
\endcode
 */
#endif

#ifndef CCP4_MMTZLIB_INC
#define CCP4_MMTZLIB_INC

#include "umtzlib.h"

#ifdef  __cplusplus
extern "C" {
#endif


/* define types for mtz records */

typedef umtzfile* mmtzfile;

/* info about mtz column */
typedef struct mmtz_column_ {
  char label[31];
  char type[2];
} mmtz_column;

/* info about mtz dataset */
typedef struct mmtz_dataset_ {
  char dname[65];
  float wavel;
} mmtz_dataset;

/* info about mtz crystal */
typedef struct mmtz_crystal_ {
  char xname[65];
  float cell[6];
  char pname[65];
} mmtz_crystal;

/* open a file for read/write (and read headers if necessary) */
mmtzfile mmtz_open( const char* filename, const char* mode );

/* close a file for read/write (and write header if necessary) */
void mmtz_close( mmtzfile file );

/* copy headers from one file to another (for append) */
void mmtz_copy_headers( mmtzfile dest, const mmtzfile src );

/* get number of reflections */
int mmtz_num_rows( const mmtzfile file );

/* get number of columns */
int mmtz_num_cols( const mmtzfile file );

/* get number of columns */
int mmtz_num_datasets( const mmtzfile file );

/* get number of symops */
int mmtz_get_num_symops( const mmtzfile file );

/* get n'th symops */
char* mmtz_get_symop( const mmtzfile file, const int isym, char* symop );

/* get base cell */
float* mmtz_get_cell( const mmtzfile file, float* cell );

/* get dataset and crystal info */
void mmtz_get_setxtl( const mmtzfile file, const int icol, mmtz_dataset* set, mmtz_crystal* xtl );

/* get column: Get the n'th column from the file, with dataset and crystal. */
int mmtz_get_column( const mmtzfile file, const int icol, mmtz_column* col, mmtz_dataset* set, mmtz_crystal* xtl );

/* add column: supply column, dataset, crystal. */
int mmtz_add_column( mmtzfile file, const mmtz_column* col, const mmtz_dataset* set, const mmtz_crystal* xtl );

/* write the initial headers to a new mtz */
void mmtz_init_headers( mmtzfile file, const char* title, const float cell[6] );

/* write sort order header to a new mtz */
void mmtz_add_sort_header( mmtzfile file, const char* sortorder );

/* write spacegroup header to a new mtz file */
void mmtz_add_syminf_header( mmtzfile file, const int nsym, const int nsymp, const char cellcode, const int ccp4_symbol, const char* hm_symbol, const char* pg_symbol );

/* write a symop header to a new mtz file */
void mmtz_add_symop_header( mmtzfile file, const char* symop );

/* get the next row of data from the file */
void mmtz_get_row( const mmtzfile file, float* fdata, int* flags );

/* write a new row of data to the file */
void mmtz_add_row( mmtzfile file, float* fdata, const int* flags );

/* skip to the n'th row in the file */
void mmtz_seek_row( const mmtzfile file, const int n );

#ifdef  __cplusplus
}
#endif


#endif
