/*cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis  
cdis 
cdis*/
/*------------------------------------------------------------------------------

getFileNames_c (c_dirname,fileNames,numfiles,match_string,status)




This routine searches for files in the configured paths 

Author: Dan Birkenheuer 
Date of development:   5/16/95
mod 6/14/95 DB increased number of filenames to 3000
mod 4/22/2002 DB increase number of filenames to 9000
mod 3/9/04 DB increased number of filesnames to 20000 for IHOP

DISCLAMER:  Note that I am by background a FORTRAN programmer.  Therefore,
though this routine is all "legal" C, it is probably written awkwardly in
the eyes of C programmers.   The purpose of this routine is to allow a fully
portable way of acquiring a list of filenames from a specified directory. 

Definitions:

c_dirname:
input directory path name (full or relative) ending in / or not (doesn't
matter)


fileNames:
an array of filenames in the calling routine defined in C as 
char filnames[20000][256]; 
or its FORTRAN equivalent:
character*256 filenames(20000)
the 19000 possible filenames is hardwired into this routine and therefore is
associated with the hardwire in the FORTRAN wrapper ment to go with this
routine (one level above) and to be the FORTRAN interface to the rest of the
FORTRAN world.  character*(*)  type dimensions are handled by this interface.


numfiles:
the number of valid files counted by this module.  Note this module will use
logic to determine possible sub directories and NOT include them in the
returned array of filenames.  Numfiles can be anything upon entry to this
routine.  It will be returned as 0 if there are no files to report or a
valid number up to 9000 counted files.

match_string:
the file matching string allowed to only contain the * wildcard and one *per
string.  This is used to filter the files identified by this routine.
defined in calling program (fortran again)
character match_string (256)

 
status: 
is the LAPS equivalent of Istatus.  1= good return, 0= failure.
currently the error condition will only be reported when an erroneous
directory path is supplied.  This routine will determine that that directory
is non-existent and report the error in this manner.  When this occurs,
numfiles will be reported as 0.




**----------------------------------------------------------------------------*/
#ifndef _getfilenames_c_
#define _getfilenames_c_
#include <config.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <unistd.h>
#include "regex.h"

#ifdef FORTRANUNDERSCORE
#define getfilenames_c getfilenames_c_
#define nstrncpy nstrncpy_
#endif
#ifdef FORTRANDOUBLEUNDERSCORE
#define getfilenames_c getfilenames_c__
#define nstrncpy nstrncpy__
#endif
#ifdef FORTRANCAPS
#define getfilenames_c GETFILENAMES_C
#define nstrncpy       NSTRNCPY
#endif


#ifdef __STDC__
void getfilenames_c ( char *c_dirname, char *fileNames, int *numfiles, 
		      char * match_string, int *status)
#else
void getfilenames_c (c_dirname,fileNames,numfiles,match_string,status )
char *c_dirname;
char *fileNames;
int *numfiles;
char *match_string;
int *status;
#endif


    {
    char match [256];
    int match_length;
    char ch [2];
    char blank[256];
    char buff [256];
    char dirname [256];
    DIR *WorkingDir;
    struct dirent *Outputofread;
    char testdirname [512];
    struct stat buf;
    int stat_status;
    int i, i2, numgood;
    regex_t preg;
    int cflags;
    char pattern [256];
    int regcomp_ret;
    regmatch_t pmatch[10];
    size_t nmatch;

    struct dirent **namelist;

#ifdef __STDC__
    int regcomp( regex_t *preg, const char *pattern, int cflags);
    int regexc ( const regex_t *preg, const char *string, size_t nmatch,
		 regmatch_t pmatch[], int eflags);
    void regfree ( regex_t *preg );
#endif


/*******************
	create blank field
*******************/

	strcpy (blank, " ");
	for (i=1; i<255; i++)
	strcat (blank, " ");

/*******************
	open the directory
*******************/  
     
	nstrncpy ( (char *) dirname, c_dirname, 255);
	nstrncpy ( (char *) match, match_string, 255);
	match_length = strlen ( match );

        *numfiles = scandir(dirname, &namelist,0,alphasort);
        if (*numfiles < 0) {
          printf("Warning from scandir - No files returned\n");
	  *numfiles=0;
	  *status = 0;
          return;
        } else {
          if( *numfiles > 20000) {
	    printf("Not all files are returned...exceeds 20000 limit.\n"); 
            printf("Num of  files = %d\n",*numfiles);
/*          free namelist entries from 20000 to end of namelist */
            for (i=20000; i < *numfiles; i++) {
              free(namelist[i]);          
            }
            *numfiles = 20000;
          }
          numgood = 0;
          for (i = 0; i < *numfiles; i++) {
	    if ( strcmp (namelist[i]->d_name,".") == 0 ) {}
	    else if ( strcmp (namelist[i]->d_name,"..") == 0 ) {}
	    else {
	      strcpy (testdirname,dirname);
	      strcat (testdirname,"/");
	      strcat(testdirname,namelist[i]->d_name);
	      stat_status = stat (testdirname,&buf);
	      if (S_ISREG( buf.st_mode ) ) {
	        strcpy (fileNames+(numgood*256) , namelist[i]->d_name);
		numgood++;
      	      } 
	    }
            free(namelist[i]);
          }
        }
        free(namelist);
	*numfiles = numgood;

/*  Commented to ....
	WorkingDir = opendir( (char *) dirname);

	if ( WorkingDir == (DIR *)NULL )
	{
	*numfiles=0;
	*status = 0;
	return;
	}

*******************
	read  the directory's files, ignore subdirectories
*******************


        Outputofread = readdir (WorkingDir);
 
	for (i=0; i<20000 && ( Outputofread != (struct dirent *)NULL ) ;)
	{
	  if(i >= 19900 )
	    printf ("WARNING numfiles nearing limit in getfilenames_c.c");
	  if ( strcmp (Outputofread->d_name,".") == 0 )
	    Outputofread = readdir (WorkingDir);
	  else if (  strcmp (Outputofread->d_name,"..") == 0 )
	    Outputofread = readdir (WorkingDir);
	
	  else
	    {
	      strcpy (testdirname,dirname);
	      strcat (testdirname,"/");
	      strcat(testdirname,Outputofread->d_name);
	      stat_status = stat (testdirname,&buf);
	      if (S_ISREG( buf.st_mode ) )
		{
		  strcpy (fileNames+(i*256) , Outputofread->d_name);
		  i++;
		  Outputofread = readdir (WorkingDir);
		}
	      else
		Outputofread = readdir (WorkingDir);
	    }
	}

	closedir (WorkingDir);
 ... here  */

	

/*******************
	recast the matching string into the regular expression matching pattn
*******************/


	strcpy (ch, " ");

	for (i = 0; i < match_length; i++)
	{
	strncpy (ch, &match[i], 1);
		if( strcmp (ch,"*") == 0 )
		{
			if (i == 0)
			strcpy(pattern,".*");
			else 
			strcat(pattern,".*");
		}


		else if( strcmp (ch,".") == 0 )
		{
			if (i == 0)
			strcpy(pattern,"\\.");
			else
			strcat(pattern,"\\.");
		}


		else if( strcmp (ch,"?") == 0 )
		{
			if (i == 0)
			strcpy(pattern,".");
			else
			strcat(pattern,".");
		}
		else
		{
			if (i == 0)
			strcpy(pattern,ch);
			else
			strcat(pattern,ch);
		}



	}

	strcat (pattern, "$");


/*******************
	set up call to regex
*******************/			

	preg.translate = 0;
	preg.fastmap = 0;
	preg.buffer = 0;
	preg.allocated = 0;
	nmatch = 0;

/*******************
	set up the pattern string given input filter
*******************/

	regcomp_ret = regcomp(&preg, pattern, REG_EXTENDED);

/*******************
	make the compare with each filename
*******************/

	i2 = 0;
	for (i=0; i< numgood; i++ )
	{
	regcomp_ret = regexec(&preg, fileNames+i*256, nmatch, pmatch, REG_EXTENDED);
	if (regcomp_ret == 0)
	{
		strcpy (buff, fileNames+i*256);
		strcpy (fileNames+i2*256, blank);
		strcpy (fileNames+i2*256 ,buff);
		i2++;
	}

	}

/*******************
	very important.... free allocated memory
*******************/

	regfree (&preg);


	*numfiles = i2;
	*status = 1;

    /*--------------------------------------------------------------------------
    ** done
    **------------------------------------------------------------------------*/
    return;
 
    }
 
#endif
  
