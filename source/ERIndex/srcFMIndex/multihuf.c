/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   bwt-based compression and indexing

   hufbzip.c
   compression and decompression using multiple huffman tables
   as in the bzip2 compressors.

   P. Ferragina & G. Manzini, 10 June 2000
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "common.h"

// ----------------------------------------------------------
// these external variables are required for the use of
// the macro bit_read_macro() and single_bit_read()
// -----------------------------------------------------------
extern uint32 Bit_buffer;
extern int  Bit_buffer_size;


#define BZ_RUNA 0
#define BZ_RUNB 1
#define BZ_N_GROUPS 6
#define BZ_G_SIZE   50
#define BZ_N_ITERS  4
#define BZ_LESSER_ICOST  0
#define BZ_GREATER_ICOST 15
#define BZ_MAX_ALPHA_SIZE 258
#define BZ_MAX_CODE_LEN    23


#define True 1
typedef unsigned short uint16;

//static __inline__ int decode_unary(void);

/* -------- arrays used by multihuf ------------ */
uchar huf_len[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];    // coding and decoding
int huf_code[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];   // coding
int rfreq[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];      // coding
int mtf_freq[BZ_MAX_ALPHA_SIZE];                // coding

uchar huf_minLens[BZ_N_GROUPS];   // decoding
int huf_limit[BZ_N_GROUPS][BZ_MAX_CODE_LEN];   // decoding
int huf_base[BZ_N_GROUPS][BZ_MAX_CODE_LEN];    // decoding
int huf_perm[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];  // decoding






/* *********************************************************
   decode a unary code read from Infile. the output is the
   # of zeroes we see before we see a 1
   1 --> 0
   01 --> 1
   001 --> 2
     etc.
   ********************************************************* */
static __inline__ int decode_unary(FM_INDEX *Infile)
{
  int t,i=0;

  do {
    single_bit_read_macro(Infile,t);
    if(t!=0) break;
    i++;
  } while(1);
  return i;
}


/* ********************************************************************
   this procedures reads a bucket from file decodes it and writes
   it to dest[] (which should be of the appropriate size).
   The decoding stops when an EOB is encountered or >= limit bytes
   have been decoded. that is, when >= limit chars have been written
   to dest the decompression terminates and the procedure returns.
   Note that more than limit chars can be
   decoded, so dest() should be large enough to contain the
   complete bucket. the procedure returns the number of chars written to
   dest (which can be less than limit (if a EOB is encountered))
   ******************************************************************** */
int multihuf_decompr(FM_INDEX *Infile,uchar *dest, int alpha_size, int limit)
{
  void hbCreateDecodeTables(int *limit,int *base,int *perm,uchar *length,
                           int minLen, int maxLen, int alphaSize );
  void out_of_mem(char *);
  void fatal_error(char *);
  // __inline__ int bit_read(int);    repalce by bit_read_macro
  int t, i, j, minLen, maxLen, len, nGroups;

  alpha_size+= 2;  // we temporarily use a larger alphabet

  // get number of groups
  bit_read_macro(Infile,nGroups,3);
  /*--- get the coding tables ---*/
  {
    int curr,uc;

    for (t = 0; t < nGroups; t++) {
      bit_read_macro(Infile,curr,5);
      for (i = 0; i < alpha_size; i++) {
	while (True) {
	  if (curr < 1 || curr > 20)
	    fatal_error("multihuf_decompr");
	  single_bit_read_macro(Infile,uc);
	  if (uc == 0) break;
	  single_bit_read_macro(Infile,uc);
	  if (uc == 0) curr++; else curr--;
	}
        huf_len[t][i] = curr;
      }
    }

    /*--- Create the Huffman decoding tables ---*/
    for (t = 0; t < nGroups; t++) {
      minLen = 32;
      maxLen = 0;
      for (i = 0; i < alpha_size; i++) {
	if (huf_len[t][i] > maxLen) maxLen = huf_len[t][i];
	if (huf_len[t][i] < minLen) minLen = huf_len[t][i];
      }
      hbCreateDecodeTables (
            &(huf_limit[t][0]),
            &(huf_base[t][0]),
            &(huf_perm[t][0]),
            &(huf_len[t][0]),
            minLen, maxLen, alpha_size
	    );
      huf_minLens[t] = minLen;
    }
  }

   /*------- uncompress data -------*/
  {
     int rle_sofar,run,next,rank,gSel,to_be_read;
     int zn,zj,zvec, *gLimit, *gPerm, *gBase;
     uchar pos[BZ_N_GROUPS], gMinlen=0;

     gLimit=gPerm=gBase=NULL;  // to avoid annoying compiler warnings
     for (i = 0; i < nGroups; i++) pos[i] = i;
     len = 0; rle_sofar=0;
     to_be_read=0;
     while (True) {
       if(to_be_read==0) {
	 to_be_read = BZ_G_SIZE;
	 rank=decode_unary(Infile);    // get mtf rank of new group
	 assert(rank<nGroups);
	 gSel=pos[rank];
	 for(j=rank;j>0;j--)  pos[j]=pos[j-1];
	 pos[0]=(uchar) gSel;
	 // get tables for this group
	 gMinlen = huf_minLens[gSel];
	 gLimit = &(huf_limit[gSel][0]);
	 gPerm = &(huf_perm[gSel][0]);
	 gBase = &(huf_base[gSel][0]);
       }
       to_be_read--;
       // get next huffman encoded char
       zn = gMinlen;
       // zvec = bit_read(zn);
       bit_read_macro(Infile,zvec,zn);
       while (zvec > gLimit[zn]) {
	 zn++;
	 // zj=bit_read(1);
	 single_bit_read_macro(Infile,zj);
	 zvec = (zvec << 1) | zj;
       };
       next = gPerm[zvec - gBase[zn]];
       // decode next
       assert(next<alpha_size);
       if(next==alpha_size-1) break;  // end of bucket
       if(next==BZ_RUNA) {            // 0 of a 1-2 encoding
	 run=1<<rle_sofar;
	 for(j=0;j<run;j++) dest[len++]=0;
	 rle_sofar++;
       }
       else if(next==BZ_RUNB) {       // 1 of a 1-2 encoding
         run = 2<<rle_sofar;
	 for(j=0;j<run;j++) dest[len++]=0;
	 rle_sofar++;
       }
       else {
	 dest[len++] = next-1;
	 rle_sofar=0;
       }
       if(len>=limit) return len; // only line added to stop when >= limit
     }                            // chars have been decoded
  }
  return len;
}


/* ************************************************************
   uncompress the bucket which starts at the current position
   of infile. the bucket is "len" bytes long and should
   be written in array dest[]
   ************************************************************ */
void uncompress_bucket_multihuf(FM_INDEX *Infile,uchar *dest, int len, int alpha_size)
{
  int int_log2(int);
  void init_bit_buffer(FM_INDEX *);
  //int bit_read(int);
  int uint_read(FM_INDEX *Infile);
  //int bit_read(FILE *,int);
   int bit_read(FM_INDEX *Infile,int n);
  void fatal_error(char *);
  void out_of_mem(char *);
  //int multihuf_decompr(uchar *, int, int);
  int multihuf_decompr(FM_INDEX *,uchar *, int, int);
  int k,j,i,aux_len,bits_x_char,rank,local_alpha_size,next;
  uchar mtf[256], inv_local_map[256];

  assert(Infile->Mtf_save==256);
  /* ---------- init ------------ */
  init_bit_buffer(Infile);         // this should not be necessary
  /* ---------- read local boolean map and compute inverse map ------ */
  local_alpha_size=0;
  for(k=0;k<alpha_size;k++)
    if(bit_read(Infile,1))
      inv_local_map[local_alpha_size++] = k;
  /* --------- read initial status of the mtf list ---------- */
  //  mtf_size = MIN(Mtf_save,local_alpha_size);
  bits_x_char = int_log2(local_alpha_size);   // read initial mtf list
  for(k=0;k<local_alpha_size;k++) {
    mtf[k] = bit_read(Infile,bits_x_char);
  }

  /* ------- decode multiple huffman codes ----- */
  //  aux = (uchar *) malloc(len*sizeof(uchar));
  //if(aux==NULL) out_of_mem("uncompress_bucket_multihuff");
  aux_len = multihuf_decompr(Infile,dest,local_alpha_size,len);
  if(aux_len!=len)
    fatal_error("Error in decompression! -uncompress_bucket_multihuf-\n");

  /* ------ decode *inplace* mtf_seq -------------------- */
  for(j=0; j<len;j++ ) {
    rank = dest[j];
    assert(rank<local_alpha_size);
    next=mtf[rank];              // decode mtf rank
    for(i=rank;i>0;i--)          // update mtf list
      mtf[i]=mtf[i-1];
    mtf[0]=next;                 // update mtf[0]
    dest[j]=inv_local_map[next]; // apply invamp
    assert(dest[j]<alpha_size);
  }


#if 0
  ---- old code to be removed soon ---
  for(i=j=0;j<len; ) {
    rank = aux[i++];
    assert(rank<=mtf_size);
    if(rank==0) {
      dest[j++]=mtf[0];       // rank zero: top of mtf list
    }
    else if(rank<mtf_size) {
      dest[j]=mtf[rank];     // decode mtf rank
      for(k=rank;k>0;k--)    // update mtf list
        mtf[k]=mtf[k-1];
      mtf[0]=dest[j++];      // update mtf[0] and j
    }
    else {                            // rank==mtf_size
      dest[j]=aux[i++];               // get char from file
      for(k=mtf_size-1;k>0;k--)       // update mtf
	mtf[k]=mtf[k-1];
      mtf[0]=dest[j++];     // update mtf[0] and j
    }
  }

  if(i!=aux_len) fatal_error("uncompress_bucket_multihuf");
  free(aux);
  /* ----- remap the bucket according to the inverse local map --- */
  for(j=0;j<len;j++) {
    c=dest[j];
    assert(c<local_alpha_size);
    dest[j]=inv_local_map[c];
    assert(dest[j]<alpha_size);
  }
#endif
}















