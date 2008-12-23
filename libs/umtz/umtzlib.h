/*! \file umtzlib.h \brief The micro-mtz low-level library.

This is a low-level library for accessing the MTZ file format. It is
not generally called by applications, but rather provides a basis for higher level libraries such as mmtzlib. */

#ifndef CCP4_UMTZLIB_INC
#define CCP4_UMTZLIB_INC

#include <math.h>
#include "library.h"

#ifdef  __cplusplus
extern "C" {
#endif


/* define mtz properties */

#define MTZRECLEN 80
#define MTZBATCHR (3*MTZRECLEN)
#define MTZBATINT 29
#define MTZBATFLT 156
#define MTZBATLEN ( MTZBATCHR + 4*(MTZBATINT+MTZBATFLT) )
#define MTZDATAOFF 20

typedef struct umtz_hdr_ {
  char data[MTZRECLEN];
} umtz_hdr;

typedef struct umtz_bat_ {
  char  cdata[3*MTZRECLEN];
  int   idata[MTZBATINT];
  float fdata[MTZBATFLT];
} umtz_bat;

/* define a list of header lines */

typedef struct umtz_list_ {
  int size, capacity;
  void* data;
} umtz_list;

/* now define a umtz file i/o object */

typedef struct umtzfile_ {
  char filename[200], mode[4];
  int iunit, ncol, nrow;
  float mnf;
  umtz_list headers;
  umtz_list history;
  umtz_list batches;
} umtzfile;


/* define unix-like i/o functions */

umtzfile* umtz_open( const char* filename, const char* mode );
void umtz_close( umtzfile* file );

/* accessor methods */

int umtz_num_cols( const umtzfile* file );
int umtz_num_rows( const umtzfile* file );
int umtz_num_head( const umtzfile* file );
int umtz_num_hist( const umtzfile* file );
float umtz_mnf( const umtzfile* file );
int umtz_ismnf( const umtzfile* file, float f );

umtz_hdr* umtz_first_head( const umtzfile* file );
umtz_hdr* umtz_first_hist( const umtzfile* file );
umtz_bat* umtz_first_bats( const umtzfile* file );
umtz_hdr* umtz_last_head( const umtzfile* file );
umtz_hdr* umtz_last_hist( const umtzfile* file );
umtz_bat* umtz_last_bats( const umtzfile* file );
int umtz_keymatch( const char* hdr, const char* key );
void umtz_add_head( umtzfile* file, const char* hdr );
void umtz_add_hist( umtzfile* file, const char* hdr );
void umtz_add_bats( umtzfile* file, const char* cdata, const int* idata, const float* fdata );

void umtz_seek_row( const umtzfile* file, const int n );
void umtz_get_row( const umtzfile* file, float* fdata );
void umtz_add_row( umtzfile* file, const float* fdata );

void umtz_get_cell( const umtzfile* file, const int ixtl, float* cell );
void umtz_cell_metric( const float* cell, float* metric );

int umtz_num_head_type( const umtzfile* file, char* head );

/* internal functions */

void umtz_make_rec( umtz_list* l, const int rec_len );
void umtz_copy_pad( char* d, const char* s, const int rec_len  );
void umtz_rewrite_headers_ranges( umtzfile* file );
void umtz_rewrite_headers_legacy( umtzfile* file );

#ifdef  __cplusplus
}
#endif

#endif
