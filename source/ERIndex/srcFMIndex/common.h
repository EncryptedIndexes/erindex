#ifndef FMCOMMON_H_
#define FMCOMMON_H_

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   bwt-based compression and indexing

   common.h
   P. Ferragina & G. Manzini, 10 June 2000
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef DEBUG
#define DEBUG 1   /* set DEBUG to 0 to remove assertions and extra checks */
#endif
#if !DEBUG
#define NDEBUG 1  /* do not compile assertions */
#endif
#include <assert.h>




#define SIZE_BASIC_PROLOGUE  55   /* #bytes except char prefix_count */

#define MIN(a, b) ((a)<=(b) ? (a) : (b))


/* **************************************************************
   Type and structure definitions

   (some data structures have been overestimated to reduce dynamic
    memory allocation and since they are used at construction time.)
   ************************************************************** */
typedef unsigned int uint32;
typedef unsigned char uchar;

typedef struct {
  int *occ;               // occ chars of compact alph in prev. superbuc.
  int alpha_size;         // actual size of alphabet in this superbucket
  uchar *bool_char_map;   // boolean map of chars occurring in this superbucket
  int *char_map;		  // serve a velocizzare i tempi di calcolo della funzione occ
  uchar inv_map_sb[256];
} bucket_lev1;


//Aggiunta da me per cachare interamente i bucket
typedef struct {
  int *occ;               // occ chars of superbucket chars in prev. buckets.
  int alpha_size;         // actual size of alphabet in this bucket
  uchar *bool_char_map;   // boolean map of chars occurring in this bucket
  uchar inv_map_b[256];
  uchar *unmtf_bucket;
  int **char_occs;  //aggiunta da me per velocizzare i tempi di calcolo della funzione occ
  int *char_num_occs;
} bucket_lev2;


typedef struct {
  uchar *text;          /* input text */
  int text_size;        /* size of text */
  int alpha_size;       /* actual size of (compact) alphabet in input text */
  int *sa;              /* suffix array */
  uchar *bwt;           /* BWT of input text */
  int bwt_eof_pos;      /* position of EOF within BWT */
  uchar bool_char_map[256];  /* ASCII[i] occurs iff bool_char_map[i] = 1 */
  uchar char_map[256];  /* ascii -> compact */
  int pfx_char_occ[256];/* entry i stores # of occ of chars 0..i-1 */
  bucket_lev1 *buclist_lev1;  /* array of num_bucs_lev1 superbuckets */
  int *start_lev2;      /* starting position of buckets */
  int *loc_occ;         /* locations of occurrences of the chosen_char */
  uchar chosen_char;    /* marked char */
  int skip;             /* one every "skip" occ of chosen_char is stored */
} bwi_input;

typedef struct {
  uchar *text;         /* input text */
  int text_size;       /* size of text */
  int alpha_size;      /* actual size of alphabet in input text */
  int *lf;             /* lf-mapping */
  uchar *bwt;          /* BWT of input text */
  int bwt_eof_pos;     /* position of EOF within BWT */
  bucket_lev1 **buclist_lev1;     /* array of superbuckets*/
  bucket_lev2 **buclist_lev2;     /* array of buckets*/
  uchar bool_char_map[256];
  int char_map[256];
  int inv_char_map[256];
  int bwt_occ[256];     /* entry i stores # of occ of chars 0..i-1 */
  int *start_lev2;      /* starting position of each buckets in compr file */
  int skip;             /* one occ of "chosen_char" every "skip"  */
  uchar chosen_char;    /* sampled char */
} bwi_out;


typedef struct FM_INDEX{
	FILE *file;
	uchar* File_start;  // byte where mapped file starts in memory
	uchar* File_end;    // byte where mapped file starts in memory
	uchar* File_pos;    // byte where read/write are currently positioned
	int Type_compression;    /* indicates type of compression adopted */
	int Is_dictionary;       /* a dictionary has to be processed */
	int Is_URL;              /* a dictionary of URLS has to be processed */
	int Is_huffword;         /* compression output of 7-bit huffword */
	int Bucket_size_lev1;    /* size of superbucket (multiple of 1K) */
	int Bucket_size_lev2;    /* size of bucket (multiple of 1K & divides _lev1) */
	int Mtf_save;            /* # of copied MTF element */
	int Num_bucs_lev1;       /* overall number of superbuckets */
	int Num_bucs_lev2;       /* overall number of buckets */
	int Start_prologue_info_sb; // starting byte of info superbuckets
	int Start_prologue_info_b;  // starting byte of info buckets
	int Start_prologue_occ;     // starting byte of char occurrences
	double Marked_char_freq;    // maximum frequency of the marked char
	//informations used for searching
	uchar Bool_map_sb[256]; // info alphabet for superbucket to be read
	uchar Inv_map_sb[256];         // inverse map for the current superbucket
	int Alpha_size_sb;             // actual size of alphabet in superbucket
	uchar Bool_map_b[256];  // info alphabet for the bucket to be read
	uchar Inv_map_b[256];   // inverse map for the current superbucket
	int Alpha_size_b;       // actual size of alphabet in bucket
	uchar Mtf[256]; // stores MTF-picture of bucket to be decompressed

	//caching system
	uchar **Cache_of_buckets;  // Array of pointers to uncompressed buckets
	int Num_bucs_in_cache;     // # buckets into cache
	int *LRU_queue_bucs;  // LRU ordered array of indices of uncmp. buckets
	int LRU_index;             // index of the least recently used bucket
	int Max_cached_buckets;    // Maximum number of cacheable buffers
	int Use_bwi_cache;       // !=0 if some cache is used
	double Cache_percentage;          // size of cache (% wrt uncompressed size)

	uint32 Bit_buffer;
	int  Bit_buffer_size;  /* number of unread/unwritten bits in Bit_buffer */
} FM_INDEX;



/* ================ "general" global variables ============ */
extern int Verbose;             /* produce verbose output */


//Provo a decommentarli per compilare compr_main.c
extern int Type_compression;    /* indicates type of compression adopted */
extern int Is_dictionary;       /* a dictionary has to be processed */
extern int Is_URL;              /* a dictionary of URLS has to be processed */
extern int Is_huffword;         /* compression output of 7-bit huffword */
extern int Bucket_size_lev1;    /* size of superbucket (multiple of 1K) */
extern int Bucket_size_lev2;    /* size of bucket (multiple of 1K & divides _lev1) */
extern int Mtf_save;            /* # of copied MTF element */
extern int Num_bucs_lev1;       /* overall number of superbuckets */
extern int Num_bucs_lev2;       /* overall number of buckets */
extern int Start_prologue_info_sb; // starting byte of info superbuckets
extern int Start_prologue_info_b;  // starting byte of info buckets
extern int Start_prologue_occ;     // starting byte of char occurrences
//// int Retrieved_occ;          // Number of retrieved explicit pos
extern double Marked_char_freq;    // maximum frequency of the marked char
extern char *Safile_name;      // name of possible sa file

/* ======= global variables for I/O ====== */
extern FILE *Infile;            /* input file */
extern FILE *Outfile;           /* output file */
extern int Infile_size;         /* size of input file */
extern int  Outfile_size;       /* size of the output file */

extern FILE *Infile_head;            /* header file of a dictionary */
extern FILE *Infile_dict;            /* dictionary file */


/* ======= global variables for file management ====== */
extern uchar *File_start;  // byte where mapped file starts in memory
extern uchar *File_end;    // byte where mapped file starts in memory
extern uchar *File_pos;    // byte where read/write are currently positioned
extern int Type_mem_ops;   // It may assume one of the three values below

// ---- word search strategy
#define WSUBSTRING 0 // arbitrary substrings
#define WPREFIX    1 // word prefix
#define WSUFFIX    2 // word suffix
#define WFULL      3 // full word

// ---- startegy for accessing file in bwsearch
#define EXT_MEM  1  // We operate via fseek, fopen, getc, putc
#define EXT_MMAP 2  // We map the file via mmap() function
#define IN_MEM   3  // We load the whole files in internal memory

// ---- algorithms for bucket compression ----
#define ARITH     1 // arithmetic coding
#define HIER3     2 // hierarchical 3 level
#define UNARY     3 // usary coding
#define MULTIH    4 // huffman with multiple tables

// ----  constants used for searching in bwi files (see search_main())
#define WHAT_CHAR_IS 1    // used to retrieve the char in a given pos
#define COUNT_CHAR_OCC 2  // used to count char-occs before a given pos
#define NULL_CHAR 0       // used to skip the count of char-occs


// *****************************************************************
// this is a macro identical (hopefully) to the function my_getc()
// *****************************************************************
extern int my_getc_macro_tmp;
#define my_getc_macro(c,f)   \
if(1) {                \
  switch (Type_mem_ops)\
    {                  \
    case EXT_MEM:      \
      my_getc_macro_tmp = getc(f->file);  \
      if(my_getc_macro_tmp==EOF) {  \
         fprintf(stderr,"Unexpected end of file -bit_read-\n"); \
         exit(1);      \
      }                \
      c= my_getc_macro_tmp;	    \
      break;           \
    case EXT_MMAP:     \
      c = *(f->File_pos++); \
      break;           \
    case IN_MEM:       \
      c = *(f->File_pos++); \
      break;           \
    default:           \
      fprintf(stderr,"Error in choosing memory management -- my_getc() --");\
      exit(1);         \
      break;           \
    }                  \
} else c=0  /* this is never executed */


/* *************************************************************************
   the following is a copy of the procedure bit_read() transformed to a macro
   to avoid the overhead of passing parameters etc
   ************************************************************************* */
extern uint32 t_macro;
#define bit_read_macro(Infile,dest,n)                                    \
{                                                                 \
  /* --- read groups of 8 bits until size>= n --- */              \
  while(Infile->Bit_buffer_size<n) {                                      \
    my_getc_macro(t_macro,Infile);                                \
    Infile->Bit_buffer |= (t_macro << (24-Infile->Bit_buffer_size));              \
    Infile->Bit_buffer_size += 8;                                         \
  }                                                               \
  /* ---- write n top bits in u ---- */                           \
  dest = Infile->Bit_buffer >> (32-n);                                    \
  /* ---- update buffer ---- */                                   \
  Infile->Bit_buffer <<= n;                                               \
  Infile->Bit_buffer_size -= n;                                           \
}
/* ************************************************************
   this macro reads a single bit from Bit_buffer.
   ************************************************************ */
#define single_bit_read_macro(Infile,dest)         \
{                                           \
  if(Infile->Bit_buffer_size==0) {                  \
    my_getc_macro(t_macro,Infile);          \
    dest = t_macro >>7;                     \
    Infile->Bit_buffer = t_macro << 25;             \
    Infile->Bit_buffer_size=7;                      \
  }                                         \
  else {                                    \
	Infile->Bit_buffer_size--;                      \
    dest = Infile->Bit_buffer >> 31;                \
    Infile->Bit_buffer <<= 1;                       \
  }                                         \
}


// -------------------------------------------------------
// finally, some prototypes common to many functions
// -------------------------------------------------------
void fatal_error(char *s);
void out_of_mem(char *s);
double getTime ( void );
int int_log2(int);
int int_pow2(int);
//int uint_read(void);
int bit_read24(FM_INDEX *Infile,int);//modifica 22/09/2017
int bit_read(FM_INDEX *Infile,int);
uchar my_getc(FM_INDEX *);
int my_fseek(FM_INDEX *, long, int);
int orig_my_fseek(FILE *f, long offset, int whence);
void init_bit_buffer(FM_INDEX *);

void read_basic_prologue(FM_INDEX *Infile,bwi_out *);
uchar *read_pattern_from_file(char *, int *);
int bwsearch(FM_INDEX *Infile,bwi_out *, uchar *, int, int *, int *);
int get_occ_pos(FM_INDEX* Infile,bwi_out *, int);
void single_search(FM_INDEX *Infile,char *pattern, int pattern_from_file);


int check_bwi_suffix(char *s);
int occ(FM_INDEX *,bwi_out *,int,uchar);
void get_info_sb(FM_INDEX *Infile,bwi_out *,int,int *);
uchar get_info_b(FM_INDEX *Infile,bwi_out *,uchar,int,int *,int);


void init_bwi_cache(FM_INDEX *);
void report_bwi_cache_usage(FM_INDEX *);
void disable_bwi_cache(FM_INDEX *);

//void my_open_file(char *);
FM_INDEX* my_open_file(char *);
void my_fclose(FM_INDEX *);
void orig_my_fclose(FILE *);

#endif
/* COMMON_H_ */


