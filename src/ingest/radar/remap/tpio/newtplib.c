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
/*
 *	Newtplib...
 *
 *	New version of TPIO (tape i/o library).
 */

#ifndef	LINT
static char newtplib_id[] = "@(#)newtplib.c	2.1	2/15/95";
#endif	/* LINT */

#include <stdio.h>
#include <fcntl.h>
#include <strings.h>
#include <sys/types.h>			/* tape include file needs it */
#include <sys/ioctl.h>			/* ditto */
#include <unistd.h>			/* for lseek */
#include <errno.h>
#include "newtplib.h"

#ifdef	RMT_IO_ALLOWED
/*
 *	For "rmt.h"
 *
 *	Setting "f_force_local" make everything really work.
 */
int	f_force_local = 0;
#include "rmtdefs.h"
#include "rmt.h"
#endif

/*
 *	Device types
 */

#define		NONE		-1
#define		DISK		0
#define		TAPE		1

struct mag 
{
	int tdrive;			/* tape/disk drive unit number */
	int mode;			/* 0 for read, >0 for write */
	int type;			/* -1 undefined, 0 disk, 1 tape */
	int errors;			/* number of i/o errors */
	int status;			/* -1 bad unit, 0 ok, 1 some error */
	int filenum;			/* disk file number */
	int eof;			/* 0 for data, 1 for eof */
	int multi;	/* -1 undef, 0 single file mode, 1 multi file mode */
	int force_rewind;		/* 1 to rewind at c_tpclose */
	char device[MAX_LEN];		/* device name */
} mag[MAX];

struct mag *pmag;

int init_flag = 0;

#ifdef	DEBUG
int last = -9999;
#endif	/* DEBUG */

char line[80];

/*
 *	initialize the variables
 */

c_init()
{
	int i;

	if (init_flag != 0)
		return;
	pmag = mag;
	for(i=0; i<MAX; i++)
	{
		pmag->tdrive = -1;
		pmag->type = NONE;
		pmag->errors = 0;
		pmag->status = -1;
		pmag->filenum = -1;
		pmag->eof = 0;
		pmag->multi = -1;
		pmag->force_rewind = 0;
		strcpy(pmag->device,"");
		pmag++;
	}
	init_flag++;
}

c_valid(iunit,caller)
int iunit;
char *caller;
{
	c_init();

	if (c_okdrive(iunit) == -1)
		return(-1);
	if (c_opendrive(iunit,caller) == -1)
		return(-1);
	return(0);
}

c_okdrive(iunit)
int iunit;
{
	if (iunit < 0 || iunit >= MAX)
	{
		printf("Illegal tape drive logical unit of %d\n",iunit);
		return(-1);
	}
	pmag = &mag[iunit];
	return(0);
}

c_opendrive(iunit,caller)
int iunit;
char *caller;
{
	if (pmag->tdrive == -1)
	{
		printf("%s: tape logical unit %d isn't open\n",caller,iunit);
		pmag->status = 1;
		return(-1);
	}
	pmag->status = 0;
	return(0);
}

c_tpopen(iunit,dev,mode)
int iunit,mode;
char *dev;
{
	char *ptr;
	int n, rc;

	c_init();

	if (c_okdrive(iunit) == -1)
		return(-1);

	if (pmag->tdrive != -1)			/* already open! */
	{
		printf("%s: tpopen: unit %d already in use!\n",dev,iunit);
		return(-1);
	}
	pmag->mode = mode;
	if (pmag->mode < 0)
	{
		printf("tpopen: %s:  invalid mode of %d\n",dev,pmag->mode);
		return(-1);
	}

	n = strlen(dev);

	if (n > MAX_LEN-1)
	{
		printf("tpopen:  %s:  name too long.\n",dev);
		return(-1);
	}
	else
		strncpy(pmag->device,dev,n);

/*
 *	Extract file/device name.
 */

	ptr = index(dev,':');
#ifndef	RMT_IO_ALLOWED
	if ( ptr )
	{
		printf("tpopen:  %s:  remote access not permitted.\n",dev);
		return(-1);
	}
#endif
	if (ptr == NULL)
		ptr = dev;
	else
		*ptr++;

/*
 *	Check for NULL file, and make it fail.  (Open would accept it.)
 */

	if (strlen(ptr) == 0)
	{
		printf("tpopen:  '%s' bad name\n",dev);
		return(-1);
	}

/*
 *	Easiest way to determine device type
 */
	if (strncmp(ptr,"/dev/",5) == 0)
	{
		pmag->type = TAPE;
		pmag->force_rewind = 0;
		if ( ptr[5] != 'n' )
		{
			pmag->force_rewind = 1;
			n -= 5;
			pmag->device[5] = 'n';
			strncpy( &pmag->device[6], &ptr[5], n );
		}
	}
#ifdef	DISK_IO_ALLOWED
	else
	{
		pmag->type = DISK;
		if (pmag->mode == 1)
			pmag->mode = O_WRONLY | O_CREAT | O_TRUNC;
	}
#else
	else
	{
		printf("tpopen:  '%s' not a tape device\n",ptr);
		return(-1);
	}
#endif

	rc = newopen(0);

	if (pmag->tdrive == -1)
	{
		if (pmag->mode == 0)
			sprintf(line,"tpopen: %s (read)",dev);
		else
			sprintf(line,"tpopen: %s (write)",dev);
		pmag->type = NONE;
		pmag->status = -1;
		pmag->eof = 0;
		pmag->multi = -1;
		pmag->force_rewind = 0;
		strcpy(pmag->device,"");
		perror(line);
		return(-1);
	}

/*
 *	I don't think this can happen, but we'll check for it anyways.
 */

	if ( rc == -1 )
	{
		printf("tpopen:  logic fault, new_open returns -1 with pmag->tdrive set to %d\n",
			pmag->tdrive);
		printf("tpopen:  recovery not possible.  Abort (with core dump).\n");
		abort();
	}

	return(0);

}

/*
 *
 *	subroutine c_tpread based on VAX version
 *
 *	arguments:
 *
 *		iunit		unit
 *		ibuf		character array
 *		ibuflen		length of ibuf
 *		lenr	> 0 	actual number of bytes read
 *			= 0	end of file or end of tape
 *			= -1	parity or some other error
 *
 */

c_tpread(iunit,ibuf,ibuflen)
int iunit,ibuflen;
char *ibuf;
{
	int lenr;
	extern errno;

	if (c_valid(iunit,"tpread") == -1)
		return(-1);

	if (pmag->eof == 1)
		newopen(1);

	pmag->status = 0;

#ifdef	RMT_IO_ALLOWED
	lenr = rmtread(pmag->tdrive,ibuf,ibuflen);
#else
	lenr = read(pmag->tdrive,ibuf,ibuflen);
#endif
	
#ifdef	DEBUG
	if (lenr != last)
	{
		printf("%s: read %d\n",pmag->device,lenr);
		last = lenr;
	}
	if (lenr == 0)				/* eof */
		printf("%s:  end of file detected\n",pmag->device);
#endif	/* DEBUG */
	if (lenr == -1)				/* some error */
	{
		sprintf(line,"tpread: %s",pmag->device);
		perror(line);
		if (errno == EINVAL)
			printf("%s: Record size of %d is too small.\n",
				pmag->device,ibuflen);
		error_count("tpread");
		pmag->status = 1;
		
	}
	if (lenr == 0)
		pmag->eof = 1;
	return(lenr);
}

/*
 *
 *	subroutine tpwrite based on VAX version
 *
 *	arguments:
 *
 *		iunit		unit
 *		ibuf		character array
 *		ibuflen		length of ibuf
 *		lenr	> 0 	actual number of bytes read
 *			= 0	end of file or end of tape
 *			= -1	parity or some other error
 *
 */

c_tpwrite(iunit,ibuf,ibuflen)
int iunit,ibuflen;
char *ibuf;
{
	int lenr;
	extern errno;

	if (c_valid(iunit,"tpwrite") == -1)
		return(-1);

	if (pmag->mode == 0)
	{
		printf("%s: write attempted on readonly drive\n",pmag->device);
		pmag->status = 1;
		return(-1);
	}

	if (pmag->eof == 1)
	{
		error_count("tpwrite");
		return(-1);
	}
	pmag->status = 0;

#ifdef	RMT_IO_ALLOWED
	lenr = rmtwrite(pmag->tdrive,ibuf,ibuflen);
#else
	lenr = write(pmag->tdrive,ibuf,ibuflen);
#endif	/* RMT_IO_ALLOWED */

#ifdef	DEBUG
	if (lenr != last)
	{
		printf("%s: write %d\n",pmag->device,lenr);
		last = lenr;
	}
	if (lenr == 0)
		printf("%s:  end of file detected\n",pmag->device);
#endif	/* DEBUG */
	if (lenr == -1)
	{
		sprintf(line,"tpwrite: %s",pmag->device);
		perror(line);
		if (errno == EINVAL)
			printf("%s: Record size of %d is too small.\n",
				pmag->device,ibuflen);
		error_count("tpwrite");
		pmag->status = 1;
	}

	if (lenr == 0)
	{
		printf("tpwrite:  %s end of medium\n",pmag->device);
		pmag->eof = 1;
	}
	return(lenr);
}

c_tpclose(iunit)
int iunit;
{

	if (c_valid(iunit,"tpclose") == -1)		/* failed */
		return(-1);
	if ( pmag->force_rewind )
		c_tprew( iunit );
#ifdef	RMT_IO_ALLOWED
	if (rmtclose(pmag->tdrive) == -1)
#else
	if (close(pmag->tdrive) == -1)
#endif
	{
		sprintf(line,"%s: tpclose",pmag->device);
		perror(line);
		pmag->status = 1;
		return(-1);
	}
	pmag->tdrive = -1;
	pmag->type = NONE;
	pmag->errors = 0;
	pmag->status = 0;
	pmag->filenum = -1;
	pmag->eof = 0;
	pmag->multi = -1;
	pmag->force_rewind = 0;
	return(0);
}

c_tprew(iunit)
int iunit;
{
	return(c_do_ioctl(iunit,"tprew",MTREW,1));
}

c_tpbsf(iunit,n)
int iunit,n;
{
	return(c_do_ioctl(iunit,"tpbsf",MTBSF,n));
}

c_tpbsr(iunit,n)
int iunit,n;
{
	return(c_do_ioctl(iunit,"tpbsr",MTBSR,n));
}

c_tpfsf(iunit,n)
int iunit,n;
{
	return(c_do_ioctl(iunit,"tpfsf",MTFSF,n));
}

c_tpfsr(iunit,n)
int iunit,n;
{
	return(c_do_ioctl(iunit,"tpfsr",MTFSR,n));
}

c_tpweof(iunit)
int iunit;
{
	int rc;
/*
 *	We now close and reopen the device so we don't trip on the 2gb limit.
 *	Closing the device forces a MTWEOF (mtio(4)).
	return(do_ioctl(iunit,"tpweof",MTWEOF,1));
 */
	if (c_valid(iunit,"tpweof") == -1)
		return(-1);
	rc = newopen(1);

	if ( rc == -1 )
	{
		printf("tpweof:  newopen(1) returned -1?  Fatal error!\n");
		abort();
	}
	return(0);
}

c_tprewoff(iunit)
int iunit;
{
	return(c_do_ioctl(iunit,"tprweoff",MTOFFL,1));
}

c_tpstatus(iunit)
int iunit;
{

	c_init();

	if (c_okdrive(iunit) == -1)
		return(-1);
	return(mag->status);
}

c_tpseek(iunit,k)
int iunit,k;
{
	int rc;

	c_init();

	if (c_okdrive(iunit) == -1)
		return(-1);
	if (pmag->type != DISK)
		return(-1);
#ifdef	RMT_IO_ALLOWED
	rc = rmtlseek(pmag->tdrive,k,1);
#else
	rc = lseek(pmag->tdrive,k,1);
#endif
	return(rc);
}

c_do_ioctl(iunit,name,k,n)
char *name;
int iunit,k,n;
{
	int rc;

	if (c_valid(iunit,name) == -1)
		return(-1);

#ifdef	DISK_IO_ALLOWED
	if ( pmag->type == DISK )
		rc = c_do_ioctl_disk( iunit, name, k, n );
	else
#endif
		rc = c_do_ioctl_tape( iunit, name, k, n );
	return( rc );
}

c_do_ioctl_tape(iunit,name,k,n)
char *name;
int iunit,k,n;
{
	struct mtop mt_command;
	int rc;

	mt_command.mt_op = k;
	mt_command.mt_count = n;
#ifdef	RMT_IO_ALLOWED
	rc = rmtioctl(pmag->tdrive,MTIOCTOP,&mt_command);
#else
	rc = ioctl(pmag->tdrive,MTIOCTOP,&mt_command);
#endif
	pmag->status = 0;
	if (rc == -1)
	{
		sprintf(line,"%s failed",name);
		perror(line);
		pmag->status = 1;
	}
	return(rc);
}

#ifdef	DISK_IO_ALLOWED
c_do_ioctl_disk(iunit,name,k,n)
char *name;
int iunit,k,n;
{
	int i, save, rc;

	pmag->status = 0;

/*
 *	Check to see if the user specified a negative number for file/record
 *	skipping.  If they did, switch the command and the number of files.
 *	This will reduce the size of the code.
 */

	if ( n < 0 )
	{
printf("command was %d n was %d\n",k,n);
		if ( k == MTBSF )
		{
			k = MTFSF;
			n = -n;
		}
		else if ( k == MTFSF )
		{
			k = MTBSF;
			n = -n;
		}
		else if ( k == MTBSR )
		{
			k = MTFSR;
			n = -n;
		}
		else if ( k == MTFSR )
		{
			k = MTBSR;
			n = -n;
		}
		else
			;				/* ignore */
printf("k now %d n %d\n",k,n);
	}

	if (k == MTREW)
	{
		pmag->filenum = -1;
		rc = newopen(1);
		if (rc == -1)
		{
			perror("tprew: failed for disk file???");
			error_count("tprew");
			pmag->status = 1;
		}
		return(rc);
	}
	else if (k == MTBSF)
	{
		pmag->filenum -= n;
		if (pmag->filenum < 0)
			pmag->filenum = 0;
/*
 *	Now, go back one more, as newopen() will add one
 */
		pmag->filenum--;
		rc = newopen(1);
		if (rc == -1)
			pmag->status = 1;
		return(rc);
	}
	else if (k == MTBSR)
	{
		printf("%s: not yet able to reposition disk files\n",name);
		pmag->status = 1;
			return(-1);
	}
	else if (k == MTFSF)
	{
		if (pmag->mode > 0)
		{
			printf("%s:  can't do while writing\n",name);
			pmag->status = 1;
			return(-1);
		}
		save = pmag->filenum;
		for(i=0; i<n; i++)
			newopen(1);
/*
 *	Now, see if we were successful, if we weren't go to end of file.
 */
		if (save + n != pmag->filenum)
		{
			printf("tpfsf:  early eot\n");
/*
 *	WARNING:  I don't think I've ever tested the rmtlseek() part of the
 *		  code.
 */
#ifdef	RMT_IO_ALLOWED
			rmtlseek(pmag->tdrive,0,2);
#else
			lseek(pmag->tdrive,0,2);
#endif
		}
		return(0);
	}
	else if (k == MTFSR)
	{
		printf("%s: not yet able to reposition disk files\n",
			name);
		pmag->status = 1;
		return(-1);
	}
	else if (k == MTWEOF)
	{
		if (pmag->mode == 0)
		{
			printf("tpweof:  readonly data.\n");
			pmag->status = 1;
			return(-1);
		}
		rc = newopen(1);
		if (rc == -1)
			pmag->status = 1;
		return(rc);
	}
	else if (k == MTOFFL)
	{
		printf("MTOFFL does not work for disk files.\n");
		pmag->status = 1;
		return(-1);
	}
	else
	{
		printf("Program error:  unknown request '%s'.\n",name);
		pmag->status = 1;
		return(-1);
	}
}
#endif

error_count(s)
char *s;
{
	pmag->errors++;
	if (pmag->errors >= MAX_ERRORS)
	{
		sprintf(line,"%s: %s:  too many i/o errors.  Exit\n",
				s,pmag->device);
		exit(1);
	}
}

newopen(mode)
int mode;
{
	int tmpfd;

/*
 *	Is this the first open?  (Called from c_tpopen().)
 */

	if (mode == 0)
	{
#ifdef	DISK_IO_ALLOWED
		if (pmag->type == DISK )
		{
			tmpfd = multi_open(1);
			if (tmpfd != -1)
				return(tmpfd);
		}
#endif

		tmpfd = multi_open(0);
		return(tmpfd);
	}

/*
 *	Later open
 */
	else
	{
/*
 *	If we are a tape, call the routine so we can close/reopen things.
 */
		if ( pmag->type == TAPE )
		{
			tmpfd = multi_open( 0 );
			return( tmpfd );
		}
/*
 *	If we aren't in multi-mode, just return so normal EOF processing occurs.
 */
		if (pmag->multi == 0)
			return(-1);
		tmpfd = multi_open(1);
		return(tmpfd);
	}
}

multi_open(n)
int n;
{
	int tmpfd;
	char filename[MAX_LEN+4];

	if (n == 0)
		strcpy(filename,pmag->device);
	else
	{
		pmag->filenum++;
		sprintf(filename,"%s.%-3.3d",pmag->device,pmag->filenum);
	}

#ifdef	DISK_IO_ALLOWED
	if ( pmag->type == DISK )
	{
#ifdef	RMT_IO_ALLOWED
		tmpfd = rmtopen(filename,pmag->mode,0644);
#else
		tmpfd = open(filename,pmag->mode,0644);
#endif
		if (tmpfd != -1)
		{
			if (pmag->tdrive != -1)
#ifdef	RMT_IO_ALLOWED
				rmtclose(pmag->tdrive);
#else
				close(pmag->tdrive);
#endif
			else
				pmag->errors = 0;
			pmag->tdrive = tmpfd;
			pmag->eof = 0;
			pmag->multi = n;
		}	
/*
 *	Undo the pmag->filenum++, if necessary
 */
		else
		{
			if (n != 0)
				pmag->filenum--;
		}
		return(tmpfd);
	}
	else
#endif
	{
		if ( pmag->tdrive != -1 )
#ifdef	RMT_IO_ALLOWED
			rmtclose( pmag->tdrive );
#else
			close( pmag->tdrive );
#endif
#ifdef	RMT_IO_ALLOWED
		pmag->tdrive = rmtopen( filename,pmag->mode,0644 );
#else
		pmag->tdrive = open( filename,pmag->mode,0644 );
#endif
		pmag->errors = 0;
		pmag->eof = 0;
		pmag->multi = n;
	}
	return( pmag->tdrive );
}
