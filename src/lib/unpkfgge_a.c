/*---------------------------------------------------------------------------*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*---------------------------------------------------------------------------*/
/*
** File :
**   unpkfgge.c
** Purpose :
**   Decode CWB NWP FGGE IIIB format grid record.
** Usage :
**   unpkfgge FGGE_binary_file
** Compile :
**   for ULTRIX : cc -DSWAPBYTE -c unpkfgge_a.c
**   for ALPHA	: cc -DSWAPBYTE -DLONG64 -c unpkfgge_a.c
**   for VMS    : cc /def=(SWAPBYTE) unpkfgge_a.c
**   for HP/UX	: cc -Aa	-c unpkfgge_a.c
**   for Linux  : gcc -DSWAPBYTE -c unpkfgge_a.c
** Restriction :
**   An internal working buffer is used for data points <= 65535
**   (Since the FGGE format use 16bits to encode the data points)
**   for HP/UX	: cc -Aa	-c unpkfgge_a.c
** Note :
**   This version is write for FORTRAN API.
** History :
**   Mar-02-1994  C.P.DZEN  Origional
**   May-30-1994  C.P.DZEN  Add trim_blank for some FORTRAN compiler
**                          can't pass a blank trimed file name.
**   Jun-06-1994  C.P.DZEN  Fix negetive header->A on 64 bits machine.
**   Feb-19-1997  C.P.DZEN  Fix record length(B value) > 65535(16bits) problem.
*/
/*---------------------------------------------------------------------------*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*---------------------------------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <string.h>
/*---------------------------------------------------------------------------*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*---------------------------------------------------------------------------*/
#if defined (FORTRANUNDERSCORE)
#   define open_fgge_file open_fgge_file_ 
#   define read_fgge_record read_fgge_record_
#   define decode_fgge_header decode_fgge_header_
#   define decode_fgge_data decode_fgge_data_
#   define close_fgge_file close_fgge_file_
#   define print_fgge_header print_fgge_header_
#   define print_fgge_data print_fgge_data_
#elif defined (FORTRANDOUBLEUNDERSCORE)
#   define open_fgge_file open_fgge_file__ 
#   define read_fgge_record read_fgge_record__
#   define decode_fgge_header decode_fgge_header__
#   define decode_fgge_data decode_fgge_data__
#   define close_fgge_file close_fgge_file__
#   define print_fgge_header print_fgge_header__
#   define print_fgge_data print_fgge_data__
#endif
/*---------------------------------------------------------------------------*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*---------------------------------------------------------------------------*/
void swap32(char *b1,int size)
{ int	loc;
  char	tmp;

  if ( (size%4) != 0) size += 4 - (size%4) ;
  if ( (size%4) != 0)
  { (void) fprintf(stderr, "Error in swap32(): size not mult of 4 bytes\n");
    exit(1);
  }
  for (loc = 0; loc < size; loc+=4)
  {
    tmp = b1[loc+0]; b1[loc+0] = b1[loc+3]; b1[loc+3] = tmp;
    tmp = b1[loc+1]; b1[loc+1] = b1[loc+2]; b1[loc+2] = tmp;
  }
}
/*---------------------------------------------------------------------------*/
void swap64(char *b1,int size)
{ int	loc;
  char	tmp;

  if ( (size%8) != 0) size += 8 - (size%8) ;
  if ( (size%8) != 0)
  { (void) fprintf(stderr, "Error in swap64(): size not mult of 8 bytes\n");
    exit(1);
  }
  for (loc = 0; loc < size; loc+=8)
  {
    tmp = b1[loc+0]; b1[loc+0] = b1[loc+7]; b1[loc+7] = tmp;
    tmp = b1[loc+1]; b1[loc+1] = b1[loc+6]; b1[loc+6] = tmp;
    tmp = b1[loc+2]; b1[loc+2] = b1[loc+5]; b1[loc+5] = tmp;
    tmp = b1[loc+3]; b1[loc+3] = b1[loc+4]; b1[loc+4] = tmp;
  }
}
/*---------------------------------------------------------------------------*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*---------------------------------------------------------------------------*/
#ifdef	LONG64

#define SWORD		 64		 /* Word size 64 bits  */
#define MASK		 ((unsigned long) 0xffffffffffffffffl)
#define G1BYTE(p,q,b)	 ((unsigned long) p << q >> (SWORD - b))
#define SWAP_BYTE(a,b)	 swap64(a,b)

#else

#define SWORD		 32		 /* Word size 32 bits  */
#define MASK		 ((unsigned) 0xffffffff)
#define G1BYTE(p,q,b)	 ((unsigned) p << q >> (SWORD - b))
#define SWAP_BYTE(a,b)	 swap32(a,b)

#endif
/*---------------------------------------------------------------------------*/
void ggbyte(p, u, q, b)
long  *p, *u, q, b;
{ register long   j, jq = q, jb = b, lb, qb;

  if (jq < SWORD)
  { j = 0; }
  else
  { j = jq / SWORD; /* number of words offset */
    jq %= SWORD;    /* odd bits of offset     */
  }

  if ((jq + jb) > SWORD)
  { qb = SWORD - jq;
    jb -= qb;
    lb = ((~(MASK << qb)) & (*(p + j))) << jb;
    jq = 0;
    j++;	    /* increment to next word */
    *u = lb + (G1BYTE(*(p + j), jq, jb));
  }
  else
  { *u = (G1BYTE(*(p + j), jq, jb)); }
}
/*---------------------------------------------------------------------------*/
void ggbytes(p, u, q, b, s, n)
long   *p, *u, q, b, s, n;
{ register long   i = 0, jp = 0;
  long		  jq = q;
  if (n > 0)
  { while (1)
    { ggbyte(p + jp, u + i, jq, b);
      if (++i == n) break;
      jq += b + s;
      jp += jq / SWORD;
      jq %= SWORD;
    }
  }
}
/*---------------------------------------------------------------------------*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*---------------------------------------------------------------------------*/
typedef struct fgge_head
  { long Q,S1,F1,T,C1,E1,M,X,S2,F2,N,C2,E2,CD,CM,KS,K,U1;
    long NW,JJ,MM,YY,GG,R,G,J,B,Z,A,U2,SN,MN,RT,DO,U3;
  } FGGE_HEAD ;
/*---------------------------------------------------------------------------*/
void open_fgge_file(char *fgge_name, int *status);
void read_fgge_record(char *obuf, int *status);
void decode_fgge_header(char *ibuf, FGGE_HEAD *header,int *status);
void decode_fgge_data(char *ibuf, FGGE_HEAD *header, float *data, int *status);
void close_fgge_file(int *status);
float rtran_mia(long mia);
void print_fgge_header(FGGE_HEAD *header);
void print_fgge_data(FGGE_HEAD *header,float *data);
/*---------------------------------------------------------------------------*/
#define   MAX_PTS    160000
/*---------------------------------------------------------------------------*/
FILE *fgge_fptr;
long wbuf[MAX_PTS];
/*---------------------------------------------------------------------------*/
char *trim_blank(char *string)
{ static char strbuf[256];
  char *start,*end,*strptr=strbuf;
  int  n;
  for(start=string ; *start==' ' && *start!=NULL ; start++) ;
  for(end=start    ; *end  !=' ' && *end  !=NULL ; end++  ) *strptr++ = *end ;
  *strptr = NULL ;
  return( strbuf );
}
/*---------------------------------------------------------------------------*/
void open_fgge_file(char *fgge_name, int *status)
{ *status = 0 ;
  if( (fgge_fptr=fopen(trim_blank(fgge_name), "r")) == (FILE *)NULL )
  { *status = -1 ;
    fprintf(stderr,"FGGE file open error : %s\n",fgge_name); }
}
/*---------------------------------------------------------------------------*/
void close_fgge_file(int *status)
{ *status = 0 ;
  fclose(fgge_fptr);
}
/*---------------------------------------------------------------------------*/
void read_fgge_record(char *obuf, int *status)
{ int  i,npts;
  long need_len,read_len;
  *status = 0 ;
  need_len = 48 ;
  read_len = fread(&obuf[0],sizeof(char),need_len,fgge_fptr);
  if( feof(fgge_fptr) ) { *status = -1 ; return ; }
  if( read_len!=need_len ) { *status = -2 ; return ; }
#ifdef SWAPBYTE
  SWAP_BYTE(obuf,need_len);
#endif

  ggbyte(obuf,&npts,240,16);
  
  /* printf("Total %d data points to read\n",npts); */

  need_len = npts * 2 ;
  read_len = fread(&obuf[48],sizeof(char),need_len,fgge_fptr);
  if( feof(fgge_fptr) ) { *status = -1 ; return ; }
  if( read_len!=need_len ) { *status = -2 ; return ; }
#ifdef SWAPBYTE
  SWAP_BYTE(&obuf[48],need_len);
#endif
}
/*---------------------------------------------------------------------------*/
void decode_fgge_header(char *ibuf, FGGE_HEAD *header, int *status)
{
  *status = 0 ;
  ;ggbyte(ibuf,&header->Q ,  0,12);
  ;ggbyte(ibuf,&header->S1, 12,12);
  ;ggbyte(ibuf,&header->F1, 24, 8);
  ;ggbyte(ibuf,&header->T , 32, 4);
  ;ggbyte(ibuf,&header->C1, 36,20);
  ;ggbyte(ibuf,&header->E1, 56, 8);  if( header->E1 > 127 ) header->E1 -= 256 ;
  ggbyte(ibuf,&header->M , 64, 4);
  ggbyte(ibuf,&header->X , 68, 8);
  ggbyte(ibuf,&header->S2, 76,12);
  ggbyte(ibuf,&header->F2, 88, 8);
  ggbyte(ibuf,&header->N , 96, 4);
  ggbyte(ibuf,&header->C2,100,20);
  ggbyte(ibuf,&header->E2,120, 8);  if( header->E2 > 127 ) header->E2 -= 256 ;
  ggbyte(ibuf,&header->CD,128, 8);
  ggbyte(ibuf,&header->CM,136, 8);
  ggbyte(ibuf,&header->KS,144, 8);
  ggbyte(ibuf,&header->K ,152, 8);
  ggbyte(ibuf,&header->U1,160,16);
  ;ggbyte(ibuf,&header->NW,176,16);
  ;ggbyte(ibuf,&header->JJ,192, 8);
  ;ggbyte(ibuf,&header->MM,200, 8);
  ;ggbyte(ibuf,&header->YY,208, 8);
  ggbyte(ibuf,&header->GG,216, 8);
  ggbyte(ibuf,&header->R ,224, 8);
  ggbyte(ibuf,&header->G ,232, 8);
  ;ggbyte(ibuf,&header->J ,240,16);
  ;ggbyte(ibuf,&header->B ,256,16);
  ggbyte(ibuf,&header->Z ,272,16);
  ;ggbyte(ibuf,&header->A ,288,32);
#ifdef  LONG64
  if( header->A > 0x3fffffffl ) header->A -= 0x100000000 ;
#endif
  ggbyte(ibuf,&header->U2,320,16);
  ;ggbyte(ibuf,&header->SN,336,16);
  ;ggbyte(ibuf,&header->MN,352, 8);
  ;ggbyte(ibuf,&header->RT,360, 8);
  ;ggbyte(ibuf,&header->DO,368, 8);
  ggbyte(ibuf,&header->U3,376, 8);
}
/*---------------------------------------------------------------------------*/
float rtran_mia(long mia)
{ int kk, i2, i3;
  float rmia;

  kk = 0x80000000;
  i2 = mia & 0x7f000000;
  i2 >>= 24;
  i3 = mia & 0xffffff;
  rmia = (float)(pow((double)16,(double)(i2-70)) * (double)i3);
  if (( mia & kk) == kk) rmia = -rmia;
  /* printf("rmia:%f, i2:%d, i3:%d, mia:%x\n",rmia,i2,i3,mia); */
  return( rmia );
}
/*---------------------------------------------------------------------------*/
void decode_fgge_data(char *ibuf, FGGE_HEAD *header, float *data, int *status)
{ float rmia,ddosn;
  int	i;
  *status = 0 ;

  rmia	= rtran_mia( header->A );
  ddosn = pow( (double)2.0 , (double)(15 - header->SN) );
  ggbytes(ibuf,wbuf,384,16,0,header->J);
  for( i=0 ; i < header->J ; i++)
  {
    if( wbuf[i] > 0x7fff ) wbuf[i] -= 0x10000 ;
    data[i] = rmia + ((double)wbuf[i] / ddosn) ;
  }
}
/*---------------------------------------------------------------------------*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*---------------------------------------------------------------------------*/
