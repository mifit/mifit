

/*                                                                          */
/* % Copyright Daresbury Laboratory 1992--1995                              */
/* % This is a CCP4 `part (i)' file for the purposes of copyright.          */
/* % See the CCP4 distribution conditions for explanation.                  */
/*                                                                          */
/* % \documentstyle[a4wide,times,noweb,makeidx]{article}                    */
/*                                                                          */
/* % \newcommand{\ac}[1]{{\rm\normalshape\sc #1}}   % acronym               */
/*                                                                          */
/* \documentclass{article}                                                  */
/* \usepackage{a4wide,times,noweb,makeidx}                                  */
/* \newcommand{\ac}[1]{\textsc{#1}}   % acronym                             */
/* \newcommand{\meta}[1]{\mbox{$\langle$\sl #1\/$\rangle$}}                 */
/* \newcommand{\ft}{\idx{Fortran}}                                          */
/* \newcommand{\idx}[1]{#1\index{#1}}                                       */
/* \newcommand{\fixme}[1]{\index{Fixme!}[{\bf Fixme!:} #1\@.]}              */
/*                                                                          */
/* \title{C library routines}\date{$ $Date: 2001/06/25 10:22:08 $ $}        */
/* \author{This version: Dave Love, Daresbury}                              */
/*                                                                          */
/* \makeindex                                                               */
/*                                                                          */
/* \noweboptions{longchunks,smallcode}                                      */
/*                                                                          */
/* \begin{document}                                                         */
/*                                                                          */
/* \maketitle                                                               */
/*                                                                          */
/* \noindent                                                                */
/* This file contains the lowest level routines for the CCP4 Program        */
/* Suite, mainly for i/o (as required by the {\tt diskio} routines) and     */
/* bit-twiddling.  It's been partly re-engineered from a non-literate C     */
/* file and isn't properly documented yet.  \fixme{It should generate user  */
/*   documentation eventually}                                              */
/* \bigskip                                                                 */
/*                                                                          */
/* \tableofcontents                                                         */
/*                                                                          */
/*                                                                          */
/* \section{Summary}                                                        */
/*                                                                          */
/* The following routines are defined:                                      */
/* \bigskip                                                                 */
/*                                                                          */
/* \noindent                                                                */
/* \begin{tabular}{ll}                                                      */
/*                                                                          */
/* Routine and  arguments &                      Purpose \\                 */
/* \hline                                                                   */
/* [[result = ccp4_ustenv(string)]]      & set an environment variable  \\     */
/* [[iunit = ccp4_qopen(filnam,istat)]]  & open random access file using [[fopen]] \\ */
/* [[ccp4_qclose(iunit)]]                & shut random access file using [[fclose]] \\ */
/* [[qmode(iunit,mode,nmcitm)]]       & change size of item in file ops. \\ */
/* [[qread(iunit,array,nitems,ier)]]  & [[fread]] from random access file \\ */
/* [[qwrite(iunit,array,nitems)]]     & [[fwrite]] to random access file \\ */
/* [[int qrarch(iunit,ipos)]]         & set up diskio number translation \\      */
/* [[qwarch(iunit,ipos)]] & write `machine stamp' to diskio file \\         */
/* [[qseek(iunit,irec,iel,lrecl)]]    & [[fseek]] within random access file \\ */
/* [[qback(iunit,lrecl)]]             & backspace within random access file \\ */
/* [[qskip(iunit,lrecl)]]             & skip forward within random access file \\ */
/* [[cqinq(iunit,lfilnm,filnam,length)]] & inquire file status on the given stream \\ */
/* [[qlocate(iunit,locate)]]          & current position within random access file  */
/* \end{tabular}                                                            */
/*                                                                          */
/*                                                                          */
/* \section{Portability}                                                    */
/*                                                                          */
/* We aim for compatibility with K\&R\index{K&R@K\&R C} and                 */
/* \index{ANSI C@\ac{ansi} C} \ac{ansi}~C\@ as well as \idx{VAX             */
/*   C}\index{C!VAX}.  One particularly annoying                            */
/* consequence is that we can't rely on [[#elif]] preprocessor              */
/* directives.  I don't know whether anything needs to be changed for the   */
/* new \idx{DEC C}, which is apparently \ac{ansi}\dots{}                    */
/*                                                                          */
/*                                                                          */
/* \section{Code}                                                           */
/*                                                                          */
/* These are the components of the code.  The \LA{}guarded code\RA{} is     */
/* executed when we've identified the platform.                             */
/*                                                                          */
/* A literate program like this is designed for human consumption.  It      */
/* comprises `chunks' of code interspersed in commentary.  Chunks,          */
/* indicated by \LA{}\dots\RA{} in code, are macro-substituted by their     */
/* definitions (starting with \LA{}\dots\endmoddef) to produce              */
/* compilable code.  The definitions of chunks may be added to later, as    */
/* inidcated by \LA{}\dots\plusendmoddef.  Chunks are cross-referenced      */
/* by their trailing tag.                                                   */
/*                                                                          */
/* <*>=                                                                     */
/* This was a literate program.  library.nw is the original source,
   from which library.c was generated by `notangle' and from which
   printable LaTeX can be produced by `noweave' if you have those
   tools.  The noweb system is available in
   anonymous@bellcore.com:pub/norman at the time of writing. */
/* \section{Global variables}                                               */
/* \subsection{Initialised variables}                                       */
/*                                                                          */
/*  Include library.h                                                       */
#include "library.h"

static void vaxF2ieeeF (union float_uint_uchar *, int);
static void ieeeF2vaxF (union float_uint_uchar *, int);
static void convexF2ieeeF (union float_uint_uchar *, int);
static void ieeeF2convexF (union float_uint_uchar *, int);

/* <global variables>=                                                      */
static char rcsid[] = "$Id: library.c,v 1.7 2001/06/25 10:22:08 ccb Exp $";
static int initialised =  0;    /* flag to initialise data and file streams */
/*  These DISKIO file modes used to include a [[b]] since they're           */
/* binary.  This caused serious lossage when reading a (scratch) file       */
/* just after writing it ([[fread]] returned [[0]] but without any error    */
/* indication).  I've no idea why.  However, I note that the \ac{ansi} C    */
/* [[b]] doesn't have any effect in \ac{posix} anyway.  In \idx{VAX} C      */
/* [[b]] means `no conversion of carriage-control information is            */
/* attempted'.  It's not clear to me whether you want this or not.  Both    */
/* possibilites {\em seem\/} to work OK \dots\ but not in \idx{DEC~C}, at   */
/* least in \idx{OpenVMS}\@.                                                */
/* This is something to watch out for on other                              */
/* systems.\index{portability!possible problem}                             */
/*                                                                          */
/* <global variables>=                                                      */
#if defined(__DECC) && defined(VMS) || defined (_WIN32)
static char *file_attribute[] = { /* DISKIO file modes */
  "wb+",   /* 'UNKNOWN'   open as 'OLD'/'NEW' check existence */
  "wb+",   /* 'SCRATCH'   open as 'OLD' and delete on closing */
  "rb+",   /* 'OLD'       file MUST exist or program halts */
  "wb+",   /* 'NEW'       create (overwrite) new file */
  "rb"     /* 'READONLY'  self explanatory */
#else
static char *file_attribute[] = {
  "w+",   /* 'UNKNOWN'   open as 'OLD'/'NEW' check existence */
  "w+",   /* 'SCRATCH'   open as 'OLD' and delete on closing */
  "r+",   /* 'OLD'       file MUST exist or program halts */
  "w+",   /* 'NEW'       create (overwrite) new file */
  "r"     /* 'READONLY'  self explanatory */
#endif
};                                                                      
/* Here is a table of bytes per item for the different i/o modes            */
/* available.  Note the \idx{machine dependencies} in here.  The            */
/* \idx{assumption} is that we have a 32-bit machine and that               */
/* [[int]]$\equiv$[[INTEGER]]([[*4]]), [[short]]$\equiv$[[INTEGER*2]],      */
/* [[float]]$\equiv$[[REAL]].                                               */
/*                                                                          */
/* <global variables>=                                                      */
static int item_sizes[] = {
  (int) sizeof (char),                                          /* 0: bytes */
  (int) sizeof (short int),                      /* 1: (integer) half words */
  (int) sizeof (float),                                   /* 2: reals/words */
  (int) sizeof (int),           /* 3: `short complex' (pairs of half words).
                                   NB int rather than 2*short since must fit
                                   into fortran integer */
  (int) 2*sizeof (float),                    /* 4: complex (pairs of words) */
  (int) sizeof (int),           /* 5: not used */
  (int) sizeof (int)            /* 6: integers */
};
/* \subsection{Uninitialised variables}                                     */
/*                                                                          */
/* <global variables>=                                                      */
static FILE *file_stream[MAXFILES];                 /* Pointer to disk file */
static char file_name[MAXFILES][MAXFLEN];      /* Pointer to disk file name */
static int  file_bytes_per_item[MAXFILES];/* Pointer to disk file item size */
static int  file_is_scratch[MAXFILES];    /* Indicates if file is 'SCRATCH' */
static int  file_last_op [MAXFILES];    /* see man fopen rd/wr combinations */
static int file_mode[MAXFILES];               /* diskio mode of each stream */
/* <global variables>=                                                      */
static uint16 nativeIT = NATIVEIT; /* machine integer type */ 
static uint16 nativeFT = NATIVEFT; /* machine float type */
static int
    Iconvert[MAXFILES],         /* integer convserion needed on read*/
    Fconvert[MAXFILES];         /* real convserion needed on read*/

/* This gets the length of a \ft{} string ([[character*]]\meta{len}         */
/* variable) \meta{s} with trailing blanks removed.  \fixme{Avoid lossage   */
/*   on null/blank string}                                                  */
/*                                                                          */
size_t ccp4_flength (char *s, int len)
{
  while (s[--len] == ' ');
  return (++len);
}
/****************************************************************************/
/* %def file_last_op file_mode                                              */
/*                                                                          */
/* \section{Internal routines}                                              */
/*                                                                          */
/* This interface to [[ccperr]] avoids mixing C and \ft{} i/o, as was       */
/* originally done.\index{error reporting}                                  */
/*                                                                          */
/* <internal routines>=                                                     */
void CCP4_fatal (const char *message)
{
  //printf(" Last system error message: %s\n",strerror(errno));
  printf(" PROGRAM:%s\n",message);
  exit(1);
 }
/* This prints a non-fatal [[message]] using the Fortran i/o.               */
/*                                                                          */
/* <internal routines>=                                                     */
void qprint (const char *message)
{
  printf ("%s\n",message);
 }
/* This reports a fatal error with a given file.                            */
/*                                                                          */
/* <internal routines>=                                                     */
void file_fatal (char *message, char *file)
{
  char *buff;
  size_t l;

  l = strlen (message) + strlen (file) + 1;
  buff = (char*)malloc (l);
  if (buff == NULL)
    CCP4_fatal ("Memory allocation failed");
  buff[0] = '\0';
  strcat (buff, message);
  strcat (buff, file);
  CCP4_fatal (buff);
}
/* \subsection{Non-\ac{ieee} floating-point conversion}                     */
/*                                                                          */
/* These conversion routines are based on \idx{HDF}, but do the             */
/* conversion in-place.  They do the obvious conversion between \idx{VAX},  */
/* \ac{ieee}\index{IEEE@\ac{ieee}} and \idx{Convex} formats implied by      */
/* the routine names.                                                       */
/*                                                                          */
/* <internal routines>=                                                     */
static void vaxF2ieeeF(union float_uint_uchar buffer[], int size)
{
  union float_uint_uchar out;
  unsigned char exp;
  int i;
  
  for (i = 0; i < size; i++) {
    exp = (buffer[i].c[1] << 1) | (buffer[i].c[0] >> 7); /* extract exponent */
    if (!exp && !buffer[i].c[1])        /* zero value */
      out.c[0] = out.c[1] = out.c[2] = out.c[3] = 0;
    else if (exp > 2) {         /* normal value */
      out.c[0] = buffer[i].c[1] - (uint8)1; /* subtracts 2 from exponent */
      /* copy mantissa, LSB of exponent */
      out.c[1] = buffer[i].c[0];
      out.c[2] = buffer[i].c[3];
      out.c[3] = buffer[i].c[2];
    } else if (exp) {           /* denormalized number */
      int shft;

      out.c[0] = buffer[i].c[1] & 0x80; /* keep sign, zero exponent */
      shft = 3 - exp;
      /* shift original mant by 1 or 2 to get denormalized mant */
      /* prefix mantissa with '1'b or '01'b as appropriate */
      out.c[1] = (uint8)((buffer[i].c[0] & 0x7f) >> shft) |
        (uint8)(0x10 << exp);
      out.c[2] = (uint8)(buffer[i].c[0] << (8-shft)) |
        (uint8)(buffer[i].c[3] >> shft);
      out.c[3] = (uint8)(buffer[i].c[3] << (8-shft)) |
        (uint8)(buffer[i].c[2] >> shft);
    } else {                    /* sign=1 -> infinity or NaN */
      out.c[0] = 0xff;          /* set exp to 255 */
      /* copy mantissa */
      out.c[1] = buffer[i].c[0] | (uint8)0x80; /* LSB of exp = 1 */
      out.c[2] = buffer[i].c[3];
      out.c[3] = buffer[i].c[2];
    }
    buffer[i] = out;            /* copy back result */
  }
}
/* <internal routines>=                                                     */
static void ieeeF2vaxF(union float_uint_uchar buffer[], int size)
{
  union float_uint_uchar out;
  unsigned char exp;
  int i;

  for (i=0; i<size; i++) {
    exp = (buffer[i].c[0]<<1) | (buffer[i].c[1]>>7); /* extract exponent */
    if (exp) {                  /* non-zero exponent */
      /* copy mantissa, last bit of exponent */
      out.c[0] = buffer[i].c[1];
      out.c[2] = buffer[i].c[3];
      out.c[3] = buffer[i].c[2];
      if (exp < 254)            /* normal value */
        out.c[1] = buffer[i].c[0] + (uint8)1; /* actually adds two to exp */
      else {                    /* infinity or NaN */
        if (exp == 254)         /* unrepresentable - OFL */
          /* set mant=0 for overflow */
          out.c[0] = out.c[1] = out.c[2] = out.c[3] = 0; 
        out.c[0] &= 0x7f;       /* set last bit of exp to 0 */
        out.c[1] = 0x80;        /* sign=1 exp=0 -> OFL or NaN.  this will raise
                                   a reserved operand exception if used. */
      }
    } else if (buffer[i].c[1] & 0x60) { /* denormalized value */
      int shft;
      
      shft = (buffer[i].c[1] & 0x40) ? 1 : 2; /* shift needed to normalize */
      /* shift mantissa */
      /* note last bit of exp set to 1 implicitly */
      out.c[0] = (uint8)(buffer[i].c[1] << shft) |
        (uint8)(buffer[i].c[2] >> (8-shft));
      out.c[3] = (uint8)(buffer[i].c[2] << shft) |
        (uint8)(buffer[i].c[3] >> (8-shft));
      out.c[2] = (uint8)(buffer[i].c[3] << shft);
      out.c[1] = (uint8)(buffer[i].c[0] & 0x80); /* sign */
      if (shft==1) {            /* set exp to 2 */
        out.c[1] |= 0x01;
        out.c[0] &= 0x7f;       /* set LSB of exp to 0 */
      }
    } else                      /* zero */
      out.c[0] = out.c[1] = out.c[2] = out.c[3] = 0;
    buffer[i] = out;            /* copy back the result */
  }
}
/* The \idx{Convex} format is like the \idx{VAX} with a different byte      */
/* order.  Convex does provide                                              */
/* \ac{ieee}$\leftrightarrow$native\index{IEEE@\ac{ieee}}                   */
/* conversion routines, but we need [[convexF2ieeeF]] anyhow.               */
/*                                                                          */
/* <internal routines>=                                                     */
static void convexF2ieeeF(union float_uint_uchar buffer[], int size)
{
  union float_uint_uchar out;
  unsigned char exp;
  int i;
  
  for (i = 0; i < size; i++) {
    exp = (buffer[i].c[0]<<1) | (buffer[i].c[1]>>7); /* extract exponent */
    if (!exp && !buffer[i].c[0])        /* zero value */
      out.c[0] = out.c[1] = out.c[2] = out.c[3] = 0;
    else if (exp > 2) {         /* normal value */
      out.c[0] = buffer[i].c[0] - (uint8)1; /* subtracts 2 from exponent */
      /* copy mantissa, LSB of exponent */
      out.c[1] = buffer[i].c[1];
      out.c[2] = buffer[i].c[2];
      out.c[3] = buffer[i].c[3];
    } else if (exp) {           /* denormalized number */
      int shft;
      
      out.c[0] = buffer[i].c[0] & 0x80; /* keep sign, zero exponent */
      shft = 3 - exp;
      /* shift original mant by 1 or 2 to get denormalized mant */
      /* prefix mantissa with '1'b or '01'b as appropriate */
      out.c[1] = (uint8)((buffer[i].c[1] & 0x7f) >> shft) |
        (uint8)(0x10 << exp);
      out.c[2] = (uint8)(buffer[i].c[1] << (8-shft)) |
        (uint8)(buffer[i].c[2] >> shft);
      out.c[3] = (uint8)(buffer[i].c[2] << (8-shft)) |
        (uint8)(buffer[i].c[3] >> shft);
    } else {                    /* sign=1 -> infinity or NaN */
      out.c[0] = 0xff;          /* set exp to 255 */
      /* copy mantissa */
      out.c[1] = buffer[i].c[1] | (uint8)0x80; /* LSB of exp = 1 */
      out.c[2] = buffer[i].c[2];
      out.c[3] = buffer[i].c[3];
    }
    buffer[i] = out;            /* copy back result */
  }
}
/* <internal routines>=                                                     */
static void ieeeF2convexF(union float_uint_uchar buffer[], int size)
{
  union float_uint_uchar out;
  unsigned char exp;
  int i;

  for (i=0; i < size; i++) {
    exp = (uint8)(buffer[i].c[0] << 1) |
      (uint8)(buffer[i].c[1] >> 7); /* extract exponent */
    if (exp) {                  /* non-zero exponent */
      /* copy mantissa, last bit of exponent */
      out.c[1] = buffer[i].c[1];
      out.c[3] = buffer[i].c[3];
      out.c[2] = buffer[i].c[2];
      if (exp < 254)            /* normal value */
        out.c[0] = buffer[i].c[0] + (uint8)1; /* actually adds two to exp */
      else {                    /* infinity or NaN */
        if (exp == 254)         /* unrepresentable - OFL */
          /* set mant=0 for overflow */
          out.c[0] = out.c[1] = out.c[2] = out.c[3] = 0; 
        out.c[1] &= 0x7f;       /* set last bit of exp to 0 */
        out.c[0] = 0x80;        /* sign=1 exp=0 -> OFL or NaN.  this will raise
                                   a reserved operand exception if used. */
      }
    } else if (buffer[i].c[1] & 0x60) { /* denormalized value */
      int shft;
      
      shft = (buffer[i].c[1] & 0x40) ? 1 : 2; /* shift needed to normalize */
      /* shift mantissa */
      /* note last bit of exp set to 1 implicitly */
      out.c[1] = (uint8)(buffer[i].c[1] << shft) |
        (uint8)(buffer[i].c[2] >> (8-shft));
      out.c[2] = (uint8)(buffer[i].c[2] << shft) |
        (uint8)(buffer[i].c[3] >> (8-shft));
      out.c[3] = (uint8)(buffer[i].c[3] << shft);
      out.c[0] = (uint8)(buffer[i].c[0] & 0x80); /* sign */
      if (shft==1) {            /* set exp to 2 */
        out.c[0] |= 0x01;
        out.c[1] &= 0x7f;       /* set LSB of exp to 0 */
      }
    } else                      /* zero */
      out.c[0] = out.c[1] = out.c[2] = out.c[3] = 0;
    buffer[i] = out;            /* copy back the result */
  }
}
/* \section{Miscellaneous routines}                                         */
/* \subsection{{\tt int ccp4_ustenv(\meta{string})}}                        */
/*                                                                          */
/* This sets an environment variable \meta{var} to \meta{val}, where the    */
/* argument \meta{string}[[==']]\meta{var}[['//'='//']]\meta{val}[[']].     */
/* This is for use by the `\idx{logical name}' mechanism for specifying     */
/* file connexions.  Note that a \idx{VMS} varsion is supplied in {\tt      */
/*   vms.for} and that there is no standard way of setting and              */
/* environment variable.  In a minimal \ac{posix} system it might be        */
/* necessary to twiddle the environment strings explicitly.                 */
/*                                                                          */
/*                                                                          */
/* <miscellaneous routines>=                                                */
#if ! defined (VMS)
/* <ustenv code>=                                                           */
int ccp4_ustenv (char *str)
{
#if defined (sgi) || defined (sun) || defined (__hpux) || \
    defined(_AIX) || defined(ultrix) || defined (__OSF1__) || \
    defined (__osf__) || defined (__FreeBSD__) || defined (linux) || \
    defined (titan) || defined (_WIN32)
  /* putenv is the POSIX.1, draft 3 proposed mechanism */
      /* ESV seems to have it in the SysVile universe */
  char *param;

  if ( (param = (char *) malloc(strlen(str) +1)) == NULL) 
    CCP4_fatal("CCP4_USTENV: Memory allocation failed");
  strcpy(param,str);
  return (putenv (param));
  /* note the necessary lack of free() */
#else
  /* setenv is not POSIX, BSD might have to use `index' */
  int setenv ();
  char *param1,*param2;

  if ( (param1 = (char *) malloc(strlen(str) +1)) == NULL) 
  CCP4_fatal("CCP4_USTENV: Memory allocation failed");
  strcpy(param1,str);
  if ((param2 = (char *) strchr(param1, '=')) == NULL)
  return (-1);
  *param2++ = '\0';
  return (setenv (param1, param2, 1));
#endif
}
#endif


/* \section{Diskio routines}                                                */
/*                                                                          */
/* \subsection{{\tt int ccp4_qopen(\meta{filename}, \meta{istat})}}         */
/* Opens \meta{filename} on diskio stream \meta{iunit}.  \meta{istat}       */
/* corresponds to the open mode given to [[qopen]], from which [[copen]]    */
/* is always called ---see diskio documentation.                            */
/* Returns stream number.                                                   */
/*                                                                          */
/* <diskio routines>=                                                       */
int ccp4_qopen (const char *filename, int istat)
{
  int iunit;

  if (! initialised) {
    /* note that array element 0 is unused -- using it produced
       complaints from mtzlib about a zero stream */
    for (iunit = 1; iunit < MAXFILES; iunit++) {
      file_stream[iunit]         = NULL;
      file_name[iunit][0]        = '\0';
      file_bytes_per_item[iunit] = item_sizes[DEFMODE];  /* default item size */
      file_is_scratch[iunit]     = 0;
      file_last_op[iunit]        = IRRELEVANT_OP;
      file_mode[iunit] = DEFMODE;
    }
    initialised = 1;
  }
  for (iunit = 1; iunit < MAXFILES; iunit++) /* Find next available stream */
    if (file_stream[iunit] == NULL) break;
  if (iunit == MAXFILES) 
    return -1;                /* return no more units flag */

  (void) strcpy (file_name[iunit], filename);
  file_last_op[iunit] = IRRELEVANT_OP;
  file_bytes_per_item[iunit] = item_sizes[DEFMODE]; /* default item size */
  file_mode[iunit] = DEFMODE;
  file_is_scratch[iunit] = (istat == 2);

/* There are complications involved with the \idx{VMS} code:                */
/* \begin{itemize}                                                          */
/* \item We want to be able to read files written by the old assembler      */
/*   library\index{VAX!assembler library} which wrote fixed-length records and */
/*   you can't do arbitrary seeks in such a file format.  Fortunately, in   */
/*   VAX C the file can be opened as StreamLF\index{StreamLF files} (as     */
/*   we want for C i/o) regardless of what the file header says.  Thanks    */
/*   to Peter Keller for suggesting this.  (This should also work for       */
/*   files ftp'd from a Unix box.);                                         */
/* \item We can't [[unlink]] the open file from the directory a             */
/*   posteriori.  Instead it's opened with the [[tmd]] RMS option as the    */
/*   assembler routines did;                                                */
/* \item Following the suggestion in the VMS 6.0 release notes about        */
/*   faster stream i/o, we use open option [["mbc=16"]] to increase the     */
/*   block size.  (This is supposed to be the default value with            */
/*   \idx{DEC C}.)                                                          */
/* \item However, the VAX C syntax for this ([[fopen]] with varargs)        */
/*   might not be supported by non-DEC compilers (although {\tt             */
/*   gcc}\index{GCC} does seem to have it).                                 */
/* \end{itemize}                                                            */
/*                                                                          */
/* <diskio routines>=                                                       */
#ifdef VMS
  if (file_is_scratch[iunit])
    file_stream[iunit] = fopen (file_name[iunit], file_attribute[istat - 1],
				"mbc=16", /* bigger blocksize */
				"fop=tmd"); /* temporary, delete on close */
  else
    file_stream[iunit] = fopen (file_name[iunit], file_attribute[istat - 1],
				"mbc=16", /* bigger blocksize */
				"ctx=stm", "mrs=0", "rat=cr", "rfm=stmlf");
  if (file_stream[iunit] == NULL)
    file_fatal ("(CCP4_QOPEN: can't open ", file_name[iunit]);
#else
# ifdef _MVS
  if (file_is_scratch[iunit]) {
    if ((file_stream[iunit] = tmpfile()) == NULL) 
      file_fatal ("CCP4_QOPEN: can't open ", file_name[iunit]);}
  else {
    file_stream[iunit] = fopen (file_name[iunit], file_attribute[istat - 1]);
    if (file_stream[iunit] == NULL)
      file_fatal ("CCP4_QOPEN: can't open ", file_name[iunit]);
  }
# else
  file_stream[iunit] = fopen (file_name[iunit], file_attribute[istat - 1]);
  if (file_stream[iunit] == NULL)
    file_fatal ("CCP4_QOPEN: can't open ", file_name[iunit]);
  if (file_is_scratch[iunit] && unlink (file_name[iunit])!=0)
    file_fatal ("CCP4_QOPEN: error unlinking ", file_name[iunit]);
# endif
#endif
  if (file_stream[iunit] == NULL) 
    return -2;                /* (catcher) return open failure flag (*/

  Iconvert[iunit] = Fconvert[iunit] = 0;
/* It seems the \idx{OpenVMS} (don't know which version) can easily lose its */
/* place in files.  Try flushing the output buffer before messing around    */
/* with [[fseek]].  (Thanks to Richard Bryan.)  N.B.: assumes [[*iunit]]!   */
/*                                                                          */
/* <OpenVMS seek fudge>=                                                    */
#if defined (__alpha) && defined (vms)
(void) fflush (file_stream[iunit]);
#endif
  if (fseek (file_stream[iunit], 0L, SEEK_SET) != 0)
    file_fatal("CCP4_QOPEN: fseek failed on", file_name[iunit]);
  
  return iunit;                /* return filestream number */
}
/* \subsection{{\tt int ccp4_qrarch (\meta{iunit}, \meta{ipos})}}           */
/*                                                                          */
/* For binary files with a well-determined structure in terms of            */
/* [[float]]s and [[int]]s we may want to set up the connected stream to    */
/* do transparent reading of files written on a machine with a different    */
/* architecture.  This is currently the case for map\index{map files} and   */
/* \idx{MTZ files} and this routine is called from \idx{mtzlib} and         */
/* \idx{maplib}.                                                            */
/*                                                                          */
/* [[qrarch]] reads the \idx{machine stamp} at {\em word\/} \meta{ipos}     */
/* for the diskio file on stream \meta{iunit} and sets up the appropriate   */
/* bit-twiddling for subsequent [[qread]]s on that stream.  The             */
/* information read from the file is returned in the                        */
/* form $\mbox{[[fileFT]]}+16\mbox{[[fileIT]]}$.  If the stamp is zero      */
/* (as it would be for files written with a previous version of the         */
/* library) we assume the file is in native format and needs no             */
/* conversion in [[qread]]; in this case the return value will be zero and  */
/* the caller can issue a warning.  [[Iconvert]] and [[Fconvert]] are       */
/* used by [[qread]] to determine the type of conversion (if any) to be     */
/* applied to integers and reals.                                           */
/*                                                                          */
/* Fudge:\index{fudge} Ian Tickle reports old VAX files which have a machine */
/* stamp which is byte-flipped from the correct VAX value, although it should */
/* always have been zero as far as I can see.  To accommodate this, set the */
/* logical \idx{NATIVEMTZ} and the machine stamp won't be read for any      */
/* input files for which [[qrarch]] is called.                              */
/*                                                                          */
/* Extra feature: logical/environment variable [[CONVERT_FROM]] may be set to one */
/* of [[BEIEEE]], [[LEIEEE]], [[VAX]] or [[CONVEXNATIVE]] to avoid reading the */
/* machine stamp and assume the file is from the stipulated archictecture   */
/* for all input MTZ and map files for which [[qrarch]] is called.          */
/*                                                                          */
/* N.B.: leaves the stream positioned just after the machine stamp.         */
/*                                                                          */
/* <diskio routines>=                                                       */
int ccp4_qrarch (int iunit, int ipos)
{
  uint16 fileFT, fileIT;        /* float and integer machine types of file */
  unsigned char mtstring[4];    /* machine stamp */
  char *native = getenv ("NATIVEMTZ");
  char *foreign = getenv ("CONVERT_FROM");

  if (native != NULL) return 0; 
  if (foreign != NULL) {
    if (strcmp (foreign, "BEIEEE") == 0) {
      mtstring[0] = DFNTF_BEIEEE | (DFNTF_BEIEEE << 4);
      mtstring[1] = 1 | (DFNTI_MBO << 4); }
    else if (strcmp (foreign, "LEIEEE") == 0) {
      mtstring[0] = DFNTF_LEIEEE | (DFNTF_LEIEEE << 4);
      mtstring[1] = 1 | (DFNTI_IBO << 4); }
    else if (strcmp (foreign, "VAX") == 0) {
      mtstring[0] = DFNTF_VAX | (DFNTF_VAX << 4);
      mtstring[1] = 1 | (DFNTI_IBO << 4); }
    else if (strcmp (foreign, "CONVEXNATIVE") == 0) {
      mtstring[0] = DFNTF_CONVEXNATIVE | (DFNTF_CONVEXNATIVE << 4);
      mtstring[1] = 1 | (DFNTI_MBO << 4); }  
  } else {
/* It seems the \idx{OpenVMS} (don't know which version) can easily lose its */
/* place in files.  Try flushing the output buffer before messing around    */
/* with [[fseek]].  (Thanks to Richard Bryan.)  N.B.: assumes [[iunit]]!   */
/*                                                                          */
/* <OpenVMS seek fudge>=                                                    */
#if defined (__alpha) && defined (vms)
(void) fflush (file_stream[iunit]);
#endif
    if ((fseek (file_stream[iunit], (long int) ((ipos)*item_sizes[2]),
                SEEK_SET) != 0))
      file_fatal ("CCP4_QRARCH: seek failed on ", file_name[iunit]);
    file_last_op[iunit] = READ_OP;
    if (fread (mtstring, (size_t) sizeof(char), (size_t) 4,
               file_stream[iunit]) != 4)
      file_fatal ("CCP4_QRARCH: can't read machine stamp in ", file_name[iunit]);
  }
  fileIT = (mtstring[1]>>4) & 0x0f;
  fileFT = (mtstring[0]>>4) & 0x0f;
  /* Record the need for conversion and what the file type is: */
  if (fileFT != 0 && fileFT != nativeFT)
    Fconvert[iunit] = fileFT;  /* else assume native */
  if (fileIT != 0 && fileIT != nativeIT)
    Iconvert[iunit] = fileIT;  /* else assume native */
  return (fileFT + (16*fileIT));
}
/* \subsection{{\tt subroutine qwarch(\meta{iunit}, \meta{ipos})}}          */
/* This is the complement of [[qrarch]], writing the native machine         */
/* architecture information (`\idx{machine stamp}') to diskio stream        */
/* \meta{iunit} at {\em word\/} \meta{ipos}.  Currently called              */
/* from \idx{mtzlib} and \idx{maplib}.                                      */
/*                                                                          */
/* The machine stamp in [[mtstring]] is four nibbles in order, indicating   */
/* complex and real format (must both be the same), integer format and      */
/* character format (currently irrelevant).  The last two bytes of          */
/* [[mtstring]] are currently unused and always zero.                       */
/*                                                                          */
/* N.B.: leaves the stream positioned just after the machine stamp.         */
/*                                                                          */
/* <diskio routines>=                                                       */
void ccp4_qwarch (int iunit, int ipos)
{
  unsigned char mtstring[4];    /* machine stamp */
/* It seems the \idx{OpenVMS} (don't know which version) can easily lose its */
/* place in files.  Try flushing the output buffer before messing around    */
/* with [[fseek]].  (Thanks to Richard Bryan.)  N.B.: assumes [[iunit]]!   */
/*                                                                          */
/* <OpenVMS seek fudge>=                                                    */
#if defined (__alpha) && defined (vms)
(void) fflush (file_stream[iunit]);
#endif
  if (fseek (file_stream[iunit], (long int) ((ipos)*item_sizes[2]),
             SEEK_SET) != 0)
    file_fatal ("CCP4_QWARCH: seek failed on ", file_name[iunit]);
  /* nibbles packed by masking and ORing: */
  mtstring[0] = nativeFT | (nativeFT << 4);
  mtstring[1] = 1 | (nativeIT << 4);
  mtstring[2] = mtstring[3] = 0;
  file_last_op[iunit] = WRITE_OP;
  if (fwrite (mtstring, (size_t) sizeof(char), (size_t) 4,
             file_stream[iunit]) != 4)
    file_fatal ("CCP4_QWARCH: can't write machine stamp to ", file_name[iunit]);
}
/* \subsection{{\tt int ccp4_qclose (\meta{iunit})}}                        */
/* Closes the file open on \idx{diskio} stream \meta{iunit}.                */
/*                                                                          */
/* <diskio routines>=                                                       */
int ccp4_qclose (int iunit)
{
  if (! initialised) 
    CCP4_fatal ("CCP4_QCLOSE: qopen not yet called");
  if (file_stream[iunit] != NULL) {
    if (fclose (file_stream[iunit]) == EOF) 
      file_fatal ("CCP4_QCLOSE: failed on ", file_name[iunit]);
    file_stream[iunit] = NULL;
  }
  file_name[iunit][0] = '\0';
  return (0);
}
/* \subsection{{\tt int ccp4_qmode (\meta{iunit}, \meta{mode})}}            */
/* Changes the \idx{diskio} \idx{access mode} for stream \meta{iunit} to    */
/* \meta{mode}.  The resulting size in bytes of items for transfer is       */
/* returned as \meta{size}.                                                 */
/*                                                                          */
/* <diskio routines>=                                                       */
int ccp4_qmode (int iunit, int mode)
{
  if (! initialised) 
    CCP4_fatal ("CCP4_QMODE: qopen  not yet called");

  if (mode >= 0 && mode <= 6 && mode != 5)
    file_bytes_per_item[iunit] = item_sizes[mode];
  else
    CCP4_fatal ("CCP4_QMODE: bad mode");
  file_mode[iunit] = mode;
  return  (file_bytes_per_item[iunit]);       /* return number of bytes/item */
}
/* \subsection{{\tt int ccp4_qread(\meta{iunit}, \meta{buffer},             */
/*     \meta{nitems})}}                                                     */
/*                                                                          */
/* Reads \meta{nitems} in the current mode (set by [[qmode]]) from diskio stream */
/* \meta{iunit} previously opened by [[qopen]](/[[copen]]) and returns      */
/* \meta{result} which is [[0]] on success or [[-1]] at EOF\@.              */
/* It aborts on an i/o error.                                               */
/* Numbers written in a foreign format will be translated if necessary if   */
/* the stream is connected to an MTZ or map file.                           */
/*                                                                          */
/* <diskio routines>=                                                       */
int ccp4_qread (int iunit, uint8 *buffer, int nitems)
{
  int i, n, result;

  if (! initialised) 
    CCP4_fatal ("CCP4_QREAD: qopen  not yet called");

  if (file_last_op[iunit] == WRITE_OP) {
/* It seems the \idx{OpenVMS} (don't know which version) can easily lose its */
/* place in files.  Try flushing the output buffer before messing around    */
/* with [[fseek]].  (Thanks to Richard Bryan.)  N.B.: assumes [[iunit]]!   */
/*                                                                          */
/* <OpenVMS seek fudge>=                                                    */
#if defined (__alpha) && defined (vms)
(void) fflush (file_stream[iunit]);
#endif
    if (fseek (file_stream[iunit], 0L, SEEK_CUR) != 0) {
      /* return (-1);*/
      file_fatal ("CCP4_QREAD: seek error on file ", file_name[iunit]);
    }
  }

  file_last_op[iunit] = READ_OP;
  //errno = 0;
  result = (int) fread (buffer, (size_t) file_bytes_per_item[iunit], 
			(size_t) nitems, file_stream[iunit]);
  if (result != nitems) {
    if (feof (file_stream[iunit])) return (-1);
    else {
      file_fatal ("CCP4_QREAD: i/o error on ", file_name[iunit]);
    }
  }
  n = result;
  /* <convert numbers if necessary>=                                          */
    switch (file_mode[iunit]) {
    case BYTE:
      break;
    case INT16:
      if (Iconvert[iunit])
        /* \subsubsection{Converting integers}                                      */
        /* The only possibility at present is byte-swapping (since we only deal     */
        /* with \idx{twos complement} integers).  The test in the following         */
        /* [[if]] could be short-circuited on this assumption.                      */
        /*                                                                          */
        /* <convert [[n]] short integers in [[buffer]]>=                            */
        {
        if ((Iconvert[iunit]==DFNTI_MBO && nativeIT==DFNTI_IBO) ||
            (Iconvert[iunit]==DFNTI_IBO && nativeIT==DFNTI_MBO)) {
          char j;
          for (i=0; i < n*2; i+=2) {
            j = buffer[i];
            buffer[i] = buffer[i+1];
            buffer[i+1] = j; } }
        else
          CCP4_fatal("CCP4_QREAD: bad file integer type in conversion");
        }
      break;
    case INT32:
      if (Iconvert[iunit])
        /* <convert [[n]] long integers in [[buffer]]>=                             */
        {
        if ((Iconvert[iunit]==DFNTI_MBO && nativeIT==DFNTI_IBO) ||
            (Iconvert[iunit]==DFNTI_IBO && nativeIT==DFNTI_MBO))
          /* <byte-swap [[n]] full words in [[buffer]]>=                              */
          {
            char j;
            for (i=0; i < n*4; i+=4) {
              j = buffer[i];
              buffer[i] = buffer[i+3];
              buffer[i+3] = j;
              j = buffer[i+1];
              buffer[i+1] = buffer[i+2];
              buffer[i+2] =j; }
          }
        else
          CCP4_fatal("CCP4_QREAD: bad file integer type in conversion");
        }
      break;
    case FLOAT32:
      if (Fconvert[iunit])
        /* \subsubsection{Converting reals}                                         */
        /* There are more possibilities than for integers\dots{}  Remember we use   */
        /* two stages and a canonical form.                                         */
        /*                                                                          */
        /* <convert [[n]] reals in [[buffer]]>=                                     */
        {
        switch (Fconvert[iunit]) {     /* get to BE IEEE */
           case DFNTF_VAX :
             vaxF2ieeeF((union float_uint_uchar *) buffer, n);
             break;   
           case DFNTF_CONVEXNATIVE :
             convexF2ieeeF((union float_uint_uchar *) buffer, n);
             break;
           case DFNTF_BEIEEE :
             break;
           case DFNTF_LEIEEE :
             /* <byte-swap [[n]] full words in [[buffer]]>=                              */
             {
               char j;
               for (i=0; i < n*4; i+=4) {
                 j = buffer[i];
                 buffer[i] = buffer[i+3];
                 buffer[i+3] = j;
                 j = buffer[i+1];
                 buffer[i+1] = buffer[i+2];
                 buffer[i+2] =j; }
             }
             break;
           default :
             CCP4_fatal("CCP4_QREAD: bad file real type in conversion");
           }
        /* We've now got a guaranteed big-endian \ac{ieee} [[buffer]].  Turn it     */
        /* into the native form if necessary.  (This could be done with             */
        /* [[#ifdef]] since [[nativeFT]] is constant, but presumably the compiler   */
        /* can spot that.)                                                          */
        /*                                                                          */
        /* <convert [[n]] reals in [[buffer]]>=                                     */
        switch (nativeFT) {
          case DFNTF_BEIEEE :
            break;                      /* done enough */
          case DFNTF_LEIEEE :
            /* <byte-swap [[n]] full words in [[buffer]]>=                              */
            {
              char j;
              for (i=0; i < n*4; i+=4) {
                j = buffer[i];
                buffer[i] = buffer[i+3];
                buffer[i+3] = j;
                j = buffer[i+1];
                buffer[i+1] = buffer[i+2];
                buffer[i+2] =j; }
            }
            break;
          case DFNTF_CONVEXNATIVE :
            ieeeF2convexF((union float_uint_uchar *) buffer, n);
            break;
          case DFNTF_VAX :
            ieeeF2vaxF((union float_uint_uchar *) buffer, n);
            break;
          default :
            CCP4_fatal("CCP4_QREAD: bad native real type in conversion");
          }
        }
      break;
    case COMP32:
      if (Fconvert[iunit]) {
        n = 2*n;                  /* pairs of ints */
        /* \subsubsection{Converting integers}                                      */
        /* The only possibility at present is byte-swapping (since we only deal     */
        /* with \idx{twos complement} integers).  The test in the following         */
        /* [[if]] could be short-circuited on this assumption.                      */
        /*                                                                          */
        /* <convert [[n]] short integers in [[buffer]]>=                            */
        {
        if ((Iconvert[iunit]==DFNTI_MBO && nativeIT==DFNTI_IBO) ||
            (Iconvert[iunit]==DFNTI_IBO && nativeIT==DFNTI_MBO)) {
          char j;
          for (i=0; i < n*2; i+=2) {
            j = buffer[i];
            buffer[i] = buffer[i+1];
            buffer[i+1] = j; } }
        else
          CCP4_fatal("CCP4_QREAD: bad file integer type in conversion");
        }
      }
      break;
    case COMP64:
      if (Fconvert[iunit]) {
        n = 2*result;                  /* pairs of reals */
        /* \subsubsection{Converting reals}                                         */
        /* There are more possibilities than for integers\dots{}  Remember we use   */
        /* two stages and a canonical form.                                         */
        /*                                                                          */
        /* <convert [[n]] reals in [[buffer]]>=                                     */
        {
        switch (Fconvert[iunit]) {     /* get to BE IEEE */
           case DFNTF_VAX :
             vaxF2ieeeF((union float_uint_uchar *) buffer, n);
             break;   
           case DFNTF_CONVEXNATIVE :
             convexF2ieeeF((union float_uint_uchar *) buffer, n);
             break;
           case DFNTF_BEIEEE :
             break;
           case DFNTF_LEIEEE :
             /* <byte-swap [[n]] full words in [[buffer]]>=                              */
             {
               char j;
               for (i=0; i < n*4; i+=4) {
                 j = buffer[i];
                 buffer[i] = buffer[i+3];
                 buffer[i+3] = j;
                 j = buffer[i+1];
                 buffer[i+1] = buffer[i+2];
                 buffer[i+2] =j; }
             }
             break;
           default :
             CCP4_fatal("CCP4_QREAD: bad file real type in conversion");
           }
        /* We've now got a guaranteed big-endian \ac{ieee} [[buffer]].  Turn it     */
        /* into the native form if necessary.  (This could be done with             */
        /* [[#ifdef]] since [[nativeFT]] is constant, but presumably the compiler   */
        /* can spot that.)                                                          */
        /*                                                                          */
        /* <convert [[n]] reals in [[buffer]]>=                                     */
        switch (nativeFT) {
          case DFNTF_BEIEEE :
            break;                      /* done enough */
          case DFNTF_LEIEEE :
            /* <byte-swap [[n]] full words in [[buffer]]>=                              */
            {
              char j;
              for (i=0; i < n*4; i+=4) {
                j = buffer[i];
                buffer[i] = buffer[i+3];
                buffer[i+3] = j;
                j = buffer[i+1];
                buffer[i+1] = buffer[i+2];
                buffer[i+2] =j; }
            }
            break;
          case DFNTF_CONVEXNATIVE :
            ieeeF2convexF((union float_uint_uchar *) buffer, n);
            break;
          case DFNTF_VAX :
            ieeeF2vaxF((union float_uint_uchar *) buffer, n);
            break;
          default :
            CCP4_fatal("CCP4_QREAD: bad native real type in conversion");
          }
        }
      }
      break;
    default:
      CCP4_fatal ("CCP4_QREAD: Bad mode");
    }

    return (result);
}
/* \subsection{{\tt int ccp4_qreadc(\meta{iunit}, \meta{buffer},            */
/*                                  \meta{nchars})}}                        */
/*                                                                          */
/* Fills [[CHARACTER]] buffer in byte mode from diskio stream               */
/* \meta{iunit} previously opened by [[qopen]](/[[copen]]) and returns      */
/*  the number of items read or [[-1]] or [[0]] on failure.                 */
/* Call it with a character substring if necessary to control the number    */
/* of bytes read.                                                           */
/*                                                                          */
/* <diskio routines>=                                                       */
int ccp4_qreadc (int iunit, char *buffer, size_t nchars)
{
  size_t result;

  if (! initialised) 
    CCP4_fatal ("CCP4_QREADC: qopen not yet called");

  if (file_last_op[iunit] == WRITE_OP) {
/* It seems the \idx{OpenVMS} (don't know which version) can easily lose its */
/* place in files.  Try flushing the output buffer before messing around    */
/* with [[fseek]].  (Thanks to Richard Bryan.)  N.B.: assumes [[iunit]]!   */
/*                                                                          */
/* <OpenVMS seek fudge>=                                                    */
#if defined (__alpha) && defined (vms)
(void) fflush (file_stream[iunit]);
#endif
   if (fseek (file_stream[iunit], 0L, SEEK_CUR) != 0) {
     /* return (-1); */
     file_fatal ("CCP4_QREADC: seek error on file ", file_name[iunit]);
   }
  }

  file_last_op[iunit] = READ_OP;

  result = fread (buffer, (size_t) item_sizes[BYTE], 
	     nchars, file_stream[iunit]);
  if (result != nchars) {
    if (feof (file_stream[iunit])) 
      return (-1);
    else 
      file_fatal ("CCP4_QREAD:C i/o error on ", file_name[iunit]);
  }
  return (result);
}
/* \subsection{{\tt int ccp4_qwrite (\meta{iunit}, \meta{buffer},           */
/*     \meta{nitems})}}                                                     */
/* This write \meta{nitems} items from \meta{buffer} to [[qopen]]ed         */
/* stream \meta{iunit} using the current mode.                              */
/* Returns number of items written, or exits                                */
/*                                                                          */
/* <diskio routines>=                                                       */
int ccp4_qwrite (int iunit, uint8 *buffer, int nitems)
{
  size_t result;

  if (! initialised) 
    CCP4_fatal ("CCP4_QWRITE: qopen not yet called");
  if (file_last_op[iunit] == READ_OP) {
/* It seems the \idx{OpenVMS} (don't know which version) can easily lose its */
/* place in files.  Try flushing the output buffer before messing around    */
/* with [[fseek]].  (Thanks to Richard Bryan.)  N.B.: assumes [[iunit]]!   */
/*                                                                          */
/* <OpenVMS seek fudge>=                                                    */
#if defined (__alpha) && defined (vms)
(void) fflush (file_stream[iunit]);
#endif
    if (fseek (file_stream[iunit], 0L, SEEK_CUR) != 0)
      file_fatal ("CCP4_QWRITE: seek failed on ", file_name[iunit]);
  }

  file_last_op[iunit] = WRITE_OP;

  result = fwrite (buffer, (size_t) file_bytes_per_item[iunit],
                    (size_t) nitems, file_stream[iunit]);
/* We don't (necessarily?)\ get a useful system error message from          */
/* [[CCP4_fatal]] if the write fails (e.g.\ in \idx{Irix}), hance the hint       */
/* about disc space.                                                        */
/*                                                                          */
/* <diskio routines>=                                                       */
  if ((int)result != (int)nitems)
    file_fatal ("CCP4_QWRITE: i/o error (may be out of disc space): ",
           file_name[iunit]);
  
  return (result);
}
/* \subsection{{\tt int ccp4_qwritc (\meta{iunit}, \meta{buffer})}}          */
/*                                                                          */
/* Writes [[CHARACTER*(*)]] \meta{buffer} to [[qopen]]ed                    */
/* stream \meta{iunit} in byte mode.                                        */
/* Returns the number of characters written.  Otherwise exits.              */
/*                                                                          */
/* <diskio routines>=                                                       */
int ccp4_qwritc (int iunit, char *buffer, size_t nchars)
{
  size_t result;

  if (! initialised) 
    CCP4_fatal ("CCP4_QWRITC: qopen not yet called");
  if (file_last_op[iunit] == READ_OP) {
/* It seems the \idx{OpenVMS} (don't know which version) can easily lose its */
/* place in files.  Try flushing the output buffer before messing around     */
/* with [[fseek]].  (Thanks to Richard Bryan.)  N.B.: assumes [[iunit]]!     */
/*                                                                           */
/* <OpenVMS seek fudge>=                                                     */
#if defined (__alpha) && defined (vms)
(void) fflush (file_stream[iunit]);
#endif
    if (fseek (file_stream[iunit], 0L, SEEK_CUR) != 0)
      file_fatal ("CCP4_QWRITC: seek failed on", file_name[iunit]);
  }

  file_last_op[iunit] = WRITE_OP;

  result = fwrite (buffer, (size_t) item_sizes[BYTE],
                    nchars, file_stream[iunit]);
  if (result != nchars) 
    file_fatal("CCP4_QWRITC: i/o error (may be out of disc space): ",
                           file_name[iunit]);
  return (result);
}
/* \subsection{{\tt long ccp4_qseek (\meta{iunit}, \meta{irec},             */
/*     \meta{iel}, \meta{lrecl})}}                                          */
/* Seeks to element \meta{iel} in record \meta{irec} in diskio stream       */
/* \meta{iunit} whose record length is \meta{lrecl}.                        */
/* Note: C convention, starts at 0                                          */
/*                                                                          */
/* <diskio routines>=                                                       */
long ccp4_qseek (int iunit, int irec, int iel, int lrecl)
{
  long int position;

  if (! initialised) 
    CCP4_fatal ("CCP4_QSEEK: qopen not yet called");
  position = (long) (lrecl*irec + iel);
  position *= (long) file_bytes_per_item[iunit];
  file_last_op[iunit] = IRRELEVANT_OP;
/* It seems the \idx{OpenVMS} (don't know which version) can easily lose its */
/* place in files.  Try flushing the output buffer before messing around    */
/* with [[fseek]].  (Thanks to Richard Bryan.)  N.B.: assumes [[iunit]]!   */
/*                                                                          */
/* <OpenVMS seek fudge>=                                                    */
#if defined (__alpha) && defined (vms)
(void) fflush (file_stream[iunit]);
#endif
  if (fseek (file_stream[iunit],position,SEEK_SET) != 0)
    file_fatal ("CCP4_QSEEK failed -- maybe corrupt file: ",file_name[iunit]);
/* follow lseek return convention                                           */
  return (position);
}
/* \subsection{{\tt void ccp4_rewind (\meta{iunit})}}                       */
/* Rewinds \meta{iunit}.                                                    */
/*                                                                          */
/* <diskio routines>=                                                       */
void ccp4_qrewind (int iunit)
{
  if (! initialised) 
    CCP4_fatal ("CCP4_REWIND: qopen not yet called");
  file_last_op[iunit] = IRRELEVANT_OP;
/* It seems the \idx{OpenVMS} (don't know which version) can easily lose its */
/* place in files.  Try flushing the output buffer before messing around    */
/* with [[fseek]].  (Thanks to Richard Bryan.)  N.B.: assumes [[iunit]]!   */
/*                                                                          */
/* <OpenVMS seek fudge>=                                                    */
#if defined (__alpha) && defined (vms)
(void) fflush (file_stream[iunit]);
#endif
  if (fseek (file_stream[iunit],0L,SEEK_SET) != 0)
    file_fatal ("CCP4_REWIND failed -- maybe corrupt file: ",file_name[iunit]);
  //errno = 0;
}
/* \subsection{{\tt long ccp4_qback (\meta{iunit}, \meta{lrecl})}}         */
/* Backspaces one record, of length \meta{lrecl} on diskio stream \meta{iunit}. */
/*                                                                          */
/* <diskio routines>=                                                       */
long ccp4_qback (int iunit, int lrecl)
{
  long int position;

  if (! initialised) 
    CCP4_fatal ("CCP4_QBACK: qopen not yet called");
  position = ftell (file_stream[iunit]) 
    - (lrecl)*file_bytes_per_item[iunit];
  file_last_op[iunit] = IRRELEVANT_OP;
/* It seems the \idx{OpenVMS} (don't know which version) can easily lose its */
/* place in files.  Try flushing the output buffer before messing around    */
/* with [[fseek]].  (Thanks to Richard Bryan.)  N.B.: assumes [[iunit]]!    */
/* Return the new position                                                  */
/*                                                                          */
/* <OpenVMS seek fudge>=                                                    */
#if defined (__alpha) && defined (vms)
(void) fflush (file_stream[iunit]);
#endif
  if (fseek (file_stream[iunit], position, SEEK_SET) != 0)
    file_fatal ("CCP4_QBACK failed on ", file_name[iunit]);
  return (position);
}
/* \subsection{{\tt long ccp4_qskip (\meta{iunit}, \meta{lrecl})}}          */
/* Skip forward 1 record of length \meta{lrecl} on diskio stream \meta{iunit}. */
/* Return new position.                                                     */
/*                                                                          */
/* <diskio routines>=                                                       */
long ccp4_qskip (int iunit, int lrecl)
{
  long int position;

  if (! initialised) 
    CCP4_fatal ("CCP4_QSKIP: qopen not yet called");
  position = ftell (file_stream[iunit]) +
    (lrecl)*file_bytes_per_item[iunit];
  file_last_op[iunit] = IRRELEVANT_OP;
/* It seems the \idx{OpenVMS} (don't know which version) can easily lose its */
/* place in files.  Try flushing the output buffer before messing around    */
/* with [[fseek]].  (Thanks to Richard Bryan.)  N.B.: assumes [[iunit]]!   */
/*                                                                          */
/* <OpenVMS seek fudge>=                                                    */
#if defined (__alpha) && defined (vms)
(void) fflush (file_stream[iunit]);
#endif
  if (fseek (file_stream[iunit],position,SEEK_SET) != 0)
    file_fatal ("CCP4_QSKIP failed on ", file_name[iunit]);
  return (position);
}
/* \subsection{{\tt long qinq (\meta{istrm}, \meta{filnam})}}               */
/* Returns \meta{length} of the file (if any) open on diskio stream         */
/* \meta{istrm}.                                                            */
/*                                                                          */
/* <diskio routines>=                                                       */
long ccp4_qinq (int iunit, char *filnam)
{
  long position, length;

  if (! initialised) 
    CCP4_fatal ("CCP4_QINQ: qopen not yet called");
  length = -1;                                    /* default return value */
  if (file_stream[iunit] == NULL) { 
    /* no unit open -- try file name */
    for (iunit = 1; iunit < MAXFILES; iunit++)
      if (! strcmp (filnam, file_name[iunit])) break;
  }

  if (iunit == MAXFILES || file_stream[iunit] == NULL) 
    return (length);                              /* return (-1) */

  file_last_op[iunit] = IRRELEVANT_OP;
  (void) fflush (file_stream[iunit]); /* flush the output stream */
#if 0
    /* checking the return value reportedly causes problems in ultrix
       under unknown circumstances... */
  if (fflush (file_stream[iunit]) != 0)
    file_fatal ("CCP4_QINQ: flush failed on ", file_name[iunit]);
#endif
  position = ftell (file_stream[iunit]);   /* remember current position */
/* It seems the \idx{OpenVMS} (don't know which version) can easily lose its */
/* place in files.  Try flushing the output buffer before messing around    */
/* with [[fseek]].  (Thanks to Richard Bryan.)  N.B.: assumes [[iunit]]!   */
/*                                                                          */
/* <OpenVMS seek fudge>=                                                    */
#if defined (__alpha) && defined (vms)
  (void) fflush (file_stream[iunit]);
#endif
  (void) fseek (file_stream[iunit],0L,SEEK_END); /* seek EOF */
  length = (long) ftell (file_stream[iunit]); /* get file size */
  if (fseek (file_stream[iunit],position,SEEK_SET) != 0) /* seek position */
    file_fatal ("CCP4_QINQ: seek failed on ", file_name[iunit]);

  return (length);
}
/* \subsection{{\tt long qlocate (\meta{iunit})}}                           */
/* Returns the current position \meta{locate} in the diskio stream \meta{iunit}. */
/*                                                                          */
/* <diskio routines>=                                                       */
long ccp4_qlocate (int iunit)
{
  if (! initialised) 
    CCP4_fatal ("CCP4_QLOCATE: qopen not yet called");
  if (file_stream[iunit] == NULL)
    return (-1);

  return ((long) ( ftell (file_stream[iunit]) / file_bytes_per_item[iunit]));
}

/* \section{`Magic' numbers}                                                */
/*                                                                          */
/* When, for instance, an $F$ is unobserved in a derivative, we might       */
/* want to give it a special value---a `\idx{magic number}'---possibly in   */
/* addition to a special value of the $\sigma$, like a negative one.        */
/* Using such a number in a calculation (by mistake, through ignoring the   */
/* value of $\sigma$, say) should not allow one to get half-sensible        */
/* results as one might if this number was $-9999$ or some such.  (There    */
/* is non-enforced connexion between the $F$ and its $\sigma$ in the MTZ    */
/* file, although one could think of adding extra columns to the file       */
/* with bit-encoded flags telling whether the $F$ in a given column was     */
/* observed.)                                                               */
/*                                                                          */
/* The obvious tactic with \ac{ieee} arithmetic is to use a \idx{NaN}       */
/* value in such situations.  Things may be set up so that we either get    */
/* an exception on using it in arithmetic or it silently propagates to all  */
/* values using it and its presence is indicated by a NaN in the output.    */
/* On a \idx{VAX} architecture we can't use NaN, but there is the           */
/* possibility of using a                                                   */
/* `reserved operand'\index{reserved operand|see{Rop}}                      */
/* (`\idx{Rop}') value,                                                     */
/* which will cause an exception (by experiment: when used for              */
/* floating-point arithmetic {\em or\/} printed, but not when assigned).    */
/* The \idx{Convex} native mode is similar, except that the Rop may be      */
/* printed (in the form {\tt Rop0x}\meta{fraction part}).                   */
/*                                                                          */
/* On, say, the \idx{IBM 370 architecture}---which we don't currently       */
/* support---anything's a valid floating point number, and the best ploy    */
/* is probably to use the largest representable number as the `magic'       */
/* value.  This would stand a good chance of raising an overflow            */
/* exception if used.  Anyhow, if such bad use of an undefined value is     */
/* made in a program due to insufficient checking by the code, it should    */
/* be spotted on the \ac{ieee} systems and the bug fixed---it's not         */
/* strictly necessary that it should cause a fatal error on all             */
/* architectures.                                                           */
/*                                                                          */
/* We need to provide a means of setting the magic number and checking      */
/* whether a given value is such.  These are architecture-dependent         */
/* bit-level operations, hence their presence in the C code.                */
/*                                                                          */
/* The suite doesn't currently use these routines, but should do soon.      */
/* \subsection{Setting a value: {\tt union float_uint_uchar                 */ 
/*  ccp4_nan()}}                                                            */
/*                                                                          */
/* [[nan]] was originally a \ft{} [[real function]] returning the value     */
/* (and actually done in 2 stages) with a subroutine implementation like    */
/* this called by the \ft{} function to avoid problems under \idx{VMS}      */
/* and native \idx{Convex}.  However, the \idx{f2c} calling convention      */
/* for a function loses in that case since it assumes a [[double]] value    */
/* returned which is cast to [[float]] with a SIGFPE, sigh.                 */
/*                                                                          */
/* <magic numbers>=                                                         */
union float_uint_uchar ccp4_nan ()
/* We have a choice of \idx{NaN} values in                                  */
/* \ac{ieee}\index{IEEE@\ac{ieee}} arithmetic.                              */
/* [[0xfffa5a5a]] is the one used by the \idx{MIPS} compilers as an         */
/* undefined value.  Note the hex constant is the same for both byte sexes! */
/*                                                                          */
/* <magic numbers>=                                                         */
#if NATIVEFT == DFNTF_BEIEEE || NATIVEFT == DFNTF_LEIEEE
#  define NAN 0xfffa5a5a
#endif
/* For \idx{Convex} native mode and \idx{VAX} use a \idx{Rop} value:        */
/*                                                                          */
/* <magic numbers>=                                                         */
#if NATIVEFT == DFNTF_CONVEXNATIVE
#  define NAN 0x80000000
#endif
#if NATIVEFT == DFNTF_VAX
#  define NAN 0x00008000
#endif
#ifndef NAN
#  error "NAN isn't defined (needs NATIVEFT)"
#endif
{
  union float_uint_uchar realnum;

  realnum.i = NAN;
  return (realnum);
}
/* \subsection{Testing a value: {\tt int ccp4_isnan(\meta{real})}}          */
/*                                                                          */
/* We want a \ft{} logical function [[ccp4_isnan]] to test whether its      */
/* argument is a \idx{NaN} or \idx{Rop}.  We have to do this by writing a C */
/* [[int]]-valued procedure and testing the returned value in the \ft{}     */
/* so that we don't have to assume how it represents logical values.  The   */
/* {\tt diskio}\index{diskio} library module provides the                   */
/* trivial interface [[CCP4_ISNAN]].                                        */
/*                                                                          */
/* <magic numbers>=                                                         */
int ccp4_isnan (union float_uint_uchar *realnum)
{
    /* In the \ac{ieee} case we actually return true both for \idx{NaN}s        */
    /* and for \idx{Infinity}; in either case the exponent is all ones---the    */
    /* fraction is zero for Infinity and non-zero for NaN\@.  The canonical     */
    /* test for a NaN is that it doesn't compare equal to itself, but we        */
    /* don't want to rely on the compiler avoiding a bogus optimisation anyhow. */
    /*                                                                          */
    /* <test for magic number>=                                                 */
    switch (nativeFT) {
     case DFNTF_BEIEEE :
     case DFNTF_LEIEEE :
       return ((realnum->i & 0x7f800000) == 0x7f800000); /* exponent all 1s */
    /* \idx{VAX} and \idx{Convex} \idx{Rop} has sign $=1$ and zero exponent     */
    /* with the appropriate byte sex---bit 15 and bits 7--14 respectively in    */
    /* the appropriate half-word (counting from 0).                             */
    /*                                                                          */
    /* <test for magic number>=                                                 */
      case DFNTF_CONVEXNATIVE :
        return ((realnum->i & 0xff800000) == 0x80000000);      
      case DFNTF_VAX :
        return ((realnum->i & 0x0000ff80) == 0x00008000);
      default :
        CCP4_fatal("CCP4_ISNAN: bad nativeFT");
        return 0;                   /* avoid compiler warning */
      }
}
/* \subsection{Absent data test for {\tt mtzlib}: {\tt void ccp4_bml        */
/*     (\meta{ncols}, \meta{cols})}}                                        */
