#if defined(FORTRANNOUNDERSCORE)

#define c_pause c_pause
#define c_open_g c_open_g
#define c_close_g c_close_g
#define c_read c_read
#define c_write_g c_write_g
#define findgrib_g findgrib_g
#define cv_to_ut cv_to_ut
#define cv_fr_ut cv_fr_ut
#define julian_date julian_date
#define calandar_date calander_date

#elif defined(FORTRANUNDERSCORE)

#define c_pause c_pause_
#define c_open_g c_open_g_
#define c_close_g c_close_g_
#define c_read c_read_
#define c_write_g c_write_g_
#define findgrib_g findgrib_g_
#define cv_to_ut cv_to_ut_
#define cv_fr_ut cv_fr_ut_
#define julian_date julian_date_
#define calandar_date calander_date_
#define c_swap4 c_swap4_
#define c_view4 c_view4_

#elif defined(FORTRANDOUBLEUNDERSCORE)

#define c_pause c_pause__
#define c_open_g c_open_g__
#define c_close_g c_close_g__
#define c_read c_read__
#define c_write_g c_write_g__
#define findgrib_g findgrib_g__
#define cv_to_ut cv_to_ut__
#define cv_fr_ut cv_fr_ut__
#define julian_date julian_date__
#define calandar_date calander_date__
#define c_swap4 c_swap4__
#define c_view4 c_view4__

#elif defined(STARDENT) || defined(CRAY)

#define c_pause C_PAUSE 
#define c_open_g C_OPEN_G
#define c_close_g C_CLOSE_G
#define c_read C_READ
#define c_write_g C_WRITE_G
#define findgrib_g FINDGRIB_G
#define cv_to_ut CV_TO_UT
#define cv_fr_ut CV_FR_UT
#define julian_date JULIAN_DATE
#define calandar_date CALANDER_DATE

#endif


/*#ifdef LITTLE_END*/
 #define BDS_LEN(bds)		((int) ((bds[0]<<16)+(bds[1]<<8)+bds[2]))
 #define BMS_LEN(bms)		((bms) == NULL ? 0 : (bms[0]<<16)+(bms[1]<<8)+bms[2])
 #define GDS_LEN(gds)		((int) ((gds[0]<<16)+(gds[1]<<8)+gds[2]))
 #define PDS_LEN(pds)		((int) ((pds[0]<<16)+(pds[1]<<8)+pds[2]))

 #define PDS_HAS_GDS(pds)        ((pds[7] & 128) != 0)
 #define PDS_HAS_BMS(pds)        ((pds[7] & 64) != 0)

 #define LEN_HEADER_PDS (28+42+100)
 #define END_LEN 4
/*#endif*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifdef MAC
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif 

#include<string.h>
#include<memory.h>
#include<errno.h>
#define TRUE 1
#define FALSE 0
extern int errno;

#define MAXFILES 256

FILE *file[MAXFILES];
int   init=0;

struct s_gribstat{
 int         bds_start;
 long pds_len;
 long gds_len;
 long bms_len;
 long end_len;
 long bds_len;
} sg;

/*---+----------------------------------------------------------------*/
int c_open_g(char *filename, char *faccess)

{
  int i;
  int done;
  
   /* initialize file pointer array */
   if(!init) {
     for(i=0; i<MAXFILES; i++) {
         file[i] = NULL;
     }
     init=1;
   }

 #ifdef DEBUG 
    printf("c_open_g called...\n\n");
    printf("filename : %s(end)\n",filename);
    printf("faccess : %s(end)\n",faccess);
    printf("press key to continue...\n");
    getchar();
#endif 

    /* find first free slot */
    done=0;
    i=0;
    while (!done) {
      if (i >= MAXFILES) {
	done=1;
      }
      else if (file[i] == NULL) {
	done=1;
      }
      else {
	i++;
      }
    }
	
    /* Check to see if we're full */
    if (i >= MAXFILES) {
     perror("io_utils.c -> c_open_g(): MAXFILES exceeded"); 
     exit(1); 
    }    

    /* Open the file */
    file[i] = fopen(filename,faccess);
 
 
  return( i );
}



/*---+----------------------------------------------------------------*/

int c_close_g(int *fd)
{
  int err;

  err=fclose(file[*fd]);
  file[*fd]=NULL;
  return(err);
}



/*---+----------------------------------------------------------------*/

int c_read(int *fbyte, int *numbytes, void *a, int *fd)

{
  int retcode;



  retcode=fseek(file[*fd],*fbyte,0);


  if(retcode != 0)
   {
     printf("C_read error - filename %s \n",file);
     return(retcode);
   }

  retcode=fread((char*)a,1,*numbytes,file[*fd]);

  if(retcode != *numbytes)
   {
     return(errno);
   }
  else
   {
     return(0);
   }
}

/*---+----------------------------------------------------------------*/

int c_write_g(int *fbyte, int *numbytes, void *a, int *fd)
{
  int retcode;
  int i = 0;

  if(*numbytes == 0 ){
     perror("asked c_write_g to write ZERO bytes!!");
     printf("Trying to write ZERO bytes!!\n");
     return -1;
  }
  retcode=fseek(file[*fd],*fbyte,2);
  if(retcode != 0)
   {
     printf("C_write error - filename \n");
     return(retcode);
   }

  retcode=fwrite((char*)a,1,*numbytes,file[*fd]);

  if(retcode != *numbytes)
   {

     return(errno);
   }
  else
   {
     return(0);
   }
}




/*---+----------------------------------------------------------------*/

#define WORD_SIZE sizeof(int)			/*assumed addresses are same size as INT*/



int c_swap4(int *numbytes, void *a)
{
  int retcode = 0;
  int i = 0;

  int j    = 0;                  	/*where the data starts*/
  int last = 0; 			/*index of last full word to swap*/
  int n_end = 0;			/*number of 'lose' bytes at end of odd length buffer*/

  char *b = (char *) a; 


  if(*numbytes < 4) perror("SWAPPING LESS THAN FOUR BYTES! ");

  if((long) b % WORD_SIZE != 0 ) printf("INPUT BUFFER NOT ALIGNED ON WORD BOUNDRY!\n");

  /*perror("INPUT BUFFER NOT ALIGNED ON WORD BOUNDRY! ");*/
 

  last  = (int) (((*numbytes / ((int) WORD_SIZE))-1) * ((int) WORD_SIZE));

  
  if(last < 0 || *numbytes < 0) { perror("OVERFLOW IN IO_UTILS.C - C_SWAP4"); exit(1); }

  n_end = *numbytes % WORD_SIZE;


/*code to align on a word boundry*/
/*while((int)a % sizeof(int)  !=  0) j++;*/

   


/* swap bytes from big/little endian to little/big endian */ 

    if(*numbytes >= 4){
      while(j <= last){

        b[j] =   b[j] ^ b[j+3]; 
        b[j+3] = b[j] ^ b[j+3]; 	/*swap j and j+3*/
        b[j] =   b[j] ^ b[j+3]; 

        b[j+1] = b[j+1] ^ b[j+2]; 
        b[j+2] = b[j+1] ^ b[j+2];		/*swap j+1 and j+2*/ 
        b[j+1] = b[j+1] ^ b[j+2]; 
     
        j += 4; 
      }
    }

   switch(n_end){
     case 3:	b[j+1] = b[j+1] ^ b[j+2]; 
      		b[j+2] = b[j+1] ^ b[j+2];		/*swap j+1 and j+2*/ 
      		b[j+1] = b[j+1] ^ b[j+2]; 

/*#ifdef DEBUG*/
      		if(b[j] != 0) perror("CANNOT SWAP THIS WORD, UNKNOWN BYTE # 4! ");
		retcode = 1;
/*#endif*/
		break;

#ifdef DEBUG
     case 2:	if((b[j]+b[j+1]) != 0) perror("CANNOT SWAP THIS WORD, UNKNOWN BYTE # 3 AND 4! ");
		retcode = 1;
		break;

     case 1: 	if(b[j] != 0) perror("CANNOT SWAP THIS WORD, UNKNOWN BYTE # 2 AND 3 AND 4! "); 
		retcode = 1;
		break;
#endif

     default:	break;
   }

    return(retcode);
}/* end c_swap4() */


/*---+----------------------------------------------------------------*/

int c_view4(int *numbytes, void *a, int *iy)
{
  int retcode = 0;
  int i = 0;
  int old_numbytes = *numbytes;
  long j = 0;                         /*where the data starts*/
  char C[3]; 

  char  *b = (char *) a; 
  float *c = (float *) a; 
  int   *d = (int *) a; 

  C[0] = C[1] = C[2] = 0;

  printf("\n\n\n");
  printf("c_view4 called...\n");
  printf("a's addr : %d\n",a);
  printf("numbytes : %d\n",*numbytes);
  printf("iy : %i\n",*iy);
  
 if(*iy == 0){
  printf("viewing as CHAR...\n");

  i = 0;
/*
  while(i < *numbytes && ((char) *b == NULL)){
    i++;
    b++;
  }
*/
  while(i < *numbytes && getchar() != 'q'){
    C[1] = (char) *b;
    printf("b[%i] = %o\t@ %d\n", i, C[1] ,b);
      i++;
      b++; 
  }
 }

 if(*iy == 1){
  printf("viewing as FLOAT...\n");

  i = 0;
  while(i < *numbytes && ((float) *(c) == 0)){ 
    i++;
    c++;
  }

  while(i < *numbytes && getchar() != 'q'){
    printf("c[%i] = %f\t@ %d\n", i, (float ) *(c), c);
    if(i < *numbytes && ((float) *c != 0)){
      i++;
      c++;
      continue;
    }
    while(i < *numbytes && ((float) *c == 0)){ 
      i++;
      c++;
    }
  }
 }

 if(*iy == 2){
  printf("viewing as INT...\n");

  i = 0;
  while(i < *numbytes && ((int) *d == 0) ){
    i++;
    d++;
  }

  while(getchar() != 'q'){
    printf("d[%i] = %i\t@ %d\n", i, (int ) *d, d);
    if(i < *numbytes && ( (int) *d != 0)){
      i++;
      d++;
      continue;
    }
    while(i < *numbytes && ((int) *(d) == 0)){
      i++;
      d++;
    }
  }
 }

     /**numbytes = old_numbytes;*/
     return(retcode);
}/* end c_view4() */


/*---+----------------------------------------------------------------*/

int findgrib_g(int *fd)
{
  int found, grib_byte, buf_counter, num;
  char buf[8200], *place_ptr, *t_ptr;

  grib_byte=0;
  found=FALSE;
  buf_counter=0;

  rewind(file[*fd]);

  while (!found)
  {
    fread(buf,sizeof(char),8192,file[*fd]);
    buf[8192]=0;
/*         replace occurrences of null with a blank in buf
           because strstr will stop when it reaches a null
*/
    t_ptr = buf;
    num = 8192;
    while((t_ptr = memchr(t_ptr,0,num)) != NULL)
    {
      *t_ptr=' ';
      num = 8192 - (t_ptr-buf);
    }
    place_ptr = strstr(buf,"GRIB");
    if(place_ptr == NULL)
    {
      buf_counter+=8192;
    }
    else
    {
      found=TRUE;
      grib_byte=(place_ptr - buf) + buf_counter;
    }
  }
  return(grib_byte);
}










/*************************************************************************
 * TIME_DATE.C: Routines to do time/date convertions. FORTRAN CALLABLE
 *		Note. Use Full year, i.e. 1991, not '91. (Unix time routines
 *		often leave off the '19'.) Watch out!
 *		Jan is month 1.
 *
 *	F Hage	Oct, 1992   NCAR/RAP
 *   Converted names for Stupid CRAY Naming conventions - 
 */

/*************************************************************************
 *	CONVERT_TO_UNIX_TIME:  Take the separate time fields 
 * 	and calculate the unix time.  FORTRAN CALLABLE:
 *  CALL CONVERT_TO_UNIX_TIME(YEAR,MONTH,DAY,HOUR,MIN,SEC,UTIME)
 *	Returns the unix time in UTIME 
 */

void cv_to_ut(year,month,day,hour,min,sec,utime)
	int	*year,*month,*day,*hour,*min,*sec,*utime;
{
	long	u_day,jday,days;
	long	u_time;

	u_day = julian_date(1,1,1970);
	jday = julian_date((*day),(*month),(*year));

	days = jday - u_day;

	*utime = (days * 86400) + (*hour * 3600) + (*min * 60) + *sec;
}


/*************************************************************************
 *	CONVERT_FROM_UNIX_TIME:  Take the unix time field 
 * 	and calculate the seperate time fields.  FORTRAN CALLABLE:
 *  CALL CONVERT_FROM_UNIX_TIME(YEAR,MONTH,DAY,HOUR,MIN,SEC,UTIME)
 *	Returns the time in YEAR,MONTH,DAY,HOUR,MIN,SEC 
 */

void cv_fr_ut(year,month,day,hour,min,sec,utime)
	int	*year,*month,*day,*hour,*min,*sec,*utime;
{
	long	u_day,j_day,days;

	u_day = julian_date(1,1,1970);

	j_day = (*utime / 86400);

	calandar_date((u_day + j_day),day,month,year);

	j_day = (*utime % 86400);
	*hour = j_day / 3600;
	*min = (j_day / 60) - (*hour * 60);
	*sec = j_day % 60;
}
 
/*************************************************************************
 *	JULIAN_DATE: Calc the Julian calandar Day Number
 *	As Taken from Computer Language- Dec 1990, pg 58
 */

int julian_date(int *day_ptr,int *month_ptr,int *year_ptr)
{
        int     day,month,year;
	int	a,b;
	double	yr_corr;

        day=*day_ptr;
	month=*month_ptr;
	year=*year_ptr;

	/* correct for negative year */
	yr_corr = (year > 0? 0.0: 0.75);
	if(month <=2) {
		year--;
		month += 12;
	}
	b=0;

	/* Cope with Gregorian Calandar reform */
	if(year * 10000.0 + month * 100.0 + day >= 15821015.0) {
		a = year / 100;
		b = 2 - a + a / 4;
	}
	
	return (int) ((365.25 * year - yr_corr) + (long) (30.6001 * (month +1)) + day + 1720994 + b);
}

/*************************************************************************
 *	CALANDAR_DATE: Calc the calandar Day from the Julian date
 *	As Taken from Computer Language- Dec 1990, pg 58
 *	Sets day,month,year as return values.
 */

calandar_date(jdate,day,month,year)
	int	jdate;
	int	*day,*month,*year;
{
	long	a,b,c,d,e,z,alph;

	z = jdate +1;

	/* Gregorian reform correction */
	if (z < 2299161) { 
		a = z; 
	} else {
		alph = (long) ((z - 1867216.25) / 36524.25);
		a = z + 1 + alph - alph / 4;
	}

	b = a + 1524;
	c = (long) ((b - 122.1) / 365.25);
	d = (long) (365.25 * c);
	e = (long) ((b - d) / 30.6001);
	*day = (int) b - d - (long) (30.6001 * e);
	*month = (int) (e < 13.5)? e - 1 : e - 13;
	*year = (int) (*month > 2.5)? (c - 4716) : c - 4715;

	 return 0;
}


int c_pause(void){
  printf("Press any key to continue..."); getchar();
}
