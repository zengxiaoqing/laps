/* Module To Interface C I/O to Fortran */
#include <config.h>
#include <stdio.h>

#if(SIZEOF_SHORT==4)
#define fint4 short
#elif(SIZEOF_INT==4)
#define fint4 int
#elif(SIZEOF_LONG==4)
#define fint4 long
#endif

#ifdef FORTRANCAPS
#define read_binary_field READ_BINARY_FIELD
#define in_to_im          IN_TO_IM
/* #define i4_to_byte        I4_TO_BYTE */
#endif

#ifdef FORTRANUNDERSCORE
#define read_binary_field read_binary_field_
#define in_to_im          in_to_im_
/* #define i4_to_byte        i4_to_byte_ */
#endif

#ifdef FORTRANDOUBLEUNDERSCORE
#define read_binary_field read_binary_field__
#define in_to_im          in_to_im__
/* #define i4_to_byte        i4_to_byte__ */
#endif

#define swap2(x) ((((x)>>8)&255)|(((x)&255)<<8))


int read_binary_field(char *data, fint4 *in_size, fint4 *out_size, fint4 *nitems,
                      char *fname, fint4 *len)
{
  FILE *fp;
  int cnt;

 
  fname[*len]='\0';
  fp = fopen(fname,"rb");
  cnt = fread(data,(*in_size),(*nitems),fp);
  fclose(fp);

  if(cnt != (*nitems)){
    fprintf(stderr,"Failed to read number of expected items in read_binary_field %d %d\n",
            *nitems,cnt);
    return -1;
  }
  cnt = in_to_im(in_size,out_size,data,nitems);

  return(cnt);
  
 
}  
  
/* translates a field of integer*n to a field of integer*m given that
   m>= n
*/

int in_to_im(fint4 *insize, fint4 *outsize,char *data, fint4 *nitems)
{
  int i,n,m;
  int *j;
  short *k;  

  if(*insize == *outsize)
    return 0;
  if(*insize > *outsize)
    return -1;

  n = (*insize);
  m = (*outsize);

  if(n==sizeof(short) && m==sizeof(fint4)){
    j = (fint4 *) data;
    k = (short *) data;
    for(i=(*nitems)-1;i>=0;i--)
#ifdef WORDS_BIGENDIAN
      j[i]=(fint4) k[i];
#else
      j[i]=(fint4) swap2(k[i]);
#endif
    return 0;
  }
  if(n==sizeof(char) && m==sizeof(fint4)){
    j = (fint4 *) data;
    for(i=(*nitems)-1;i>=0;i--)
      j[i]=(fint4) data[i];
    return 0;
  }
  fprintf(stderr,"Failed to find type match in in_to_im %d %d\n",m,n);
  return -1;
}

int in_to_fm(fint4 *insize, fint4 *outsize,char *data, fint4 *nitems)
{
  int j,n,m;
  short *s;
  int *i;
  float *f;  
  double *d;

  if(*insize == *outsize)
    return 0;
  if(*insize > *outsize)
    return -1;

  n = (*insize);
  m = (*outsize);

  if(n==sizeof(short) && m==sizeof(float)){
    s = (short *) data;
    f = (float *) data;
    for(j=(*nitems)-1;j>=0;j--)
      f[j]=(float) s[j];
    return 0;
  }
  if(n==sizeof(char) && m==sizeof(float)){
    f = (float *) data;
    for(j=(*nitems)-1;j>=0;j--)
      f[j]=(float) data[j];
    return 0;
  }
  if(n==sizeof(short) && m==sizeof(double)){
    s = (short *) data;
    d = (double *) data;
    for(j=(*nitems)-1;j>=0;j--)
      d[j]=(double) s[j];
    return 0;
  }
  if(n==sizeof(char) && m==sizeof(double)){
    d = (double *) data;
    for(j=(*nitems)-1;j>=0;j--)
      d[j]=(double) data[j];
    return 0;
  }
  fprintf(stderr,"Failed to find type match in in_to_fm %d %d\n",m,n);
  return -1;
}
/*
char i4_to_byte(fint4 *i4)
{
  if((int) (*i4) < 0 || (int) (*i4) > 256){
    fprintf(stderr,"Value out of range in i4_to_byte %d\n",(int) *i4);
  }
    
  return((char) (*i4));
}
*/
