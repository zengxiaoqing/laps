/*
 *  Copyright 1993, University Corporation for Atmospheric Research
 *  See ../COPYRIGHT file for copying and redistribution conditions.
 *
 *  This version of pqinsert.c has been modified to support the -k command-line
 *  argument.
 *
 *  Glen Pankow	     2 Sep 97	    1.1	    Added the -k command-line argument;
 *					    renumbered the RCS Id back to 1.1.
 *  Jim Edwards      9 June 99              Ported to gcc added support for LAPS Makefile
 */
/* $Id$ */

/* 
 * Convert files to ldm "products" and insert in local que
 */
#include "config.h"
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <signal.h>
#ifndef NO_MMAP
#include <sys/mman.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <rpc/rpc.h>
#include <errno.h>
#include <time.h>
#include "paths.h"
#include "ldm.h"
#include "pq.h"
#include "atofeedt.h"
#include "ldmprint.h"
#include "inetutil.h"
#include "ulog.h"
#include "md5.h"

#ifdef NO_ATEXIT
#include "atexit.h"
#endif

static int verbose = 0;
	/* N.B.: assumes hostname doesn't change during program execution :-) */
static char myname[HOSTNAMESIZE];
static feedtypet feedtype = EXP;
static const char *pqfname = DEFAULT_QUEUE;
static char *key = (char *)NULL;
pqueue *pq = NULL;


static void
usage(
	char *av0 /*  id string */
)
{
	(void)fprintf(stderr,
		"Usage: %s [options] filename ...\n\tOptions:\n", av0);
	(void)fprintf(stderr,
		"\t-v		verbose, tell me about each product\n");
	(void)fprintf(stderr,
		"\t-l logfile	log to a file rather than stderr\n");
	(void)fprintf(stderr,
		"\t-q queue     default \"%s\"\n", DEFAULT_QUEUE);
	(void)fprintf(stderr,
		"\t-s seqno     set initial product sequence number to \"seqno\", defaults to 0\n");
	(void)fprintf(stderr,
		"\t-f feedtype  assert your feed type as \"feedtype\", defaults to \"EXP\"\n");
	(void)fprintf(stderr,
		"\t-k key       use \"key\" instead of \"filename\" for the key\n",
		"\t             (only one \"filename\" processed per run\n");
	exit(1);
}


void
cleanup(void)
{
	if(pq)
	{
		(void) pq_close(pq);
		pq = NULL;
	}
	(void) closeulog();
}


static void
signal_handler(
	int sig
)
{
#ifdef SVR3SIGNALS
	/* 
	 * Some systems reset handler to SIG_DFL upon entry to handler.
	 * In that case, we reregister our handler.
	 */
	(void) signal(sig, signal_handler);
#endif
    switch(sig) {
      case SIGINT :
         exit(1);
      case SIGTERM :
         exit(1);
      case SIGPIPE :
         udebug("SIGPIPE");
         exit(1);
    }
	udebug("signal_handler: unhandled signal: %d", sig);
}


static void
set_sigactions(void)
{
#ifndef NO_POSIXSIGNALS
	struct sigaction sigact;

	sigemptyset(&sigact.sa_mask);
	sigact.sa_flags = 0;

	/* Ignore these */
	sigact.sa_handler = SIG_IGN;
	(void) sigaction(SIGHUP, &sigact, NULL);
	(void) sigaction(SIGPIPE, &sigact, NULL);
	(void) sigaction(SIGALRM, &sigact, NULL);
	(void) sigaction(SIGCHLD, &sigact, NULL);

	/* Handle these */
#ifdef SA_RESTART	/* SVR4, 4.3+ BSD */
	/* usually, restart system calls */
	sigact.sa_flags |= SA_RESTART;
#endif
	sigact.sa_handler = signal_handler;
	(void) sigaction(SIGTERM, &sigact, NULL);
	(void) sigaction(SIGPIPE, &sigact, NULL);
	/* Don't restart after interrupt */
	sigact.sa_flags = 0;
#ifdef SA_INTERRUPT	/* SunOS 4.x */
	sigact.sa_flags |= SA_INTERRUPT;
#endif
	(void) sigaction(SIGINT, &sigact, NULL);
#else
	
	(void) signal(SIGHUP, SIG_IGN);
	(void) signal(SIGPIPE, SIG_IGN);
	(void) signal(SIGALRM, SIG_IGN);
	(void) signal(SIGCHLD, SIG_IGN);

	(void) signal(SIGTERM, signal_handler);
	(void) signal(SIGPIPE, signal_handler);
	(void) signal(SIGINT, signal_handler);
#endif
}


#ifdef NO_MMAP
static int
fd_md5(MD5_CTX *md5ctxp, int fd, off_t st_size, signaturet signature)
{
	int nread;
	char buf[8192];
	MD5Init(md5ctxp);

	for(; st_size > 0; st_size -= nread )
	{
		nread = read(fd, buf, sizeof(buf));
		if(nread <= 0)
		{
			serror("fd_md5: read");
			return -1;
		} /* else */
		MD5Update(md5ctxp, buf, nread);
	}

	MD5Final(signature, md5ctxp);
	return 0;
}
#else
static int
mm_md5(MD5_CTX *md5ctxp, void *vp, size_t sz, signaturet signature)
{
	MD5Init(md5ctxp);

	MD5Update(md5ctxp, vp, sz);

	MD5Final(signature, md5ctxp);
	return 0;
}
#endif


main(
	int ac,
	char *av[]
)
{
	char *progname = av[0];
	char *logfname;
	int status;
	int seq_start = 0;

	logfname = "-";

	/*
	 * Check the environment for some options.
	 * May be overridden by command line switches below.
	 */
	{
		const char *ldmpqfname = getenv("LDMPQFNAME");
		if(ldmpqfname != NULL)
			pqfname = ldmpqfname;
	}

	{
	extern int optind;
	extern int opterr;
	extern char *optarg;
	int ch;
	int logmask = (LOG_MASK(LOG_ERR) | LOG_MASK(LOG_NOTICE));

	opterr = 1;

	while ((ch = getopt(ac, av, "vxl:q:f:s:k:")) != EOF)
		switch (ch) {
		case 'v':
			logmask |= LOG_MASK(LOG_INFO);
			break;
		case 'x':
			logmask |= LOG_MASK(LOG_DEBUG);
			break;
		case 'l':
			logfname = optarg;
			break;
		case 'q':
			pqfname = optarg;
			break;
		case 's':
			seq_start = atoi(optarg);
			break;
		case 'f':
			feedtype = atofeedtypet(optarg);
			if(feedtype == NONE)
			{
			    fprintf(stderr, "Unknown feedtype \"%s\"\n", optarg);
				usage(progname);	
			}
			break;
		case 'k':
			key = optarg;
			break;
		case '?':
			usage(progname);
			break;
		}

	ac -= optind; av += optind ;

	if(ac < 1) usage(progname);
	(void) setulogmask(logmask);
	}

	/*
	 * Set up error logging
	 */
	(void) openulog(ubasename(progname), LOG_NOTIME, LOG_LDM, logfname);

	/*
	 * register exit handler
	 */
	if(atexit(cleanup) != 0)
	{
		serror("atexit");
		exit(1);
	}

	/*
	 * set up signal handlers
	 */
	set_sigactions();

	/*
 	 * who am i, anyway
	 */
	(void) strcpy(myname, ghostname());

	/*
 	 * open the product queue
	 */
	if(status = pq_open(pqfname, PQ_DEFAULT, &pq))
	{
		uerror("pq_open: \"%s\" failed: %s",
			pqfname, status > 0 ? strerror(status) :
					"Internal error");
		exit(2);
	}


	{
	char *filename;
        char *newid;
	int fd;
	struct stat statb;
	product prod;
#ifdef NO_MMAP
	pqe_index index;
#endif
	MD5_CTX *md5ctxp = NULL;

	/*
	 * Allocate an MD5 context
	 */
	md5ctxp = new_MD5_CTX();
	if(md5ctxp == NULL)
	{
		serror("new_md5_CTX failed");
		exit(6);
	}


	/* These members are constant over the loop. */
	prod.info.origin = myname;
	prod.info.feedtype = feedtype;

	for(prod.info.seqno = seq_start ; ac > 0 ;
			 av++, ac--, prod.info.seqno++)
	{
		filename = *av;

		fd = open(filename, O_RDONLY, 0);
		if(fd == -1)
		{
			serror("open: %s", filename);
			continue;
		}

		if( fstat(fd, &statb) == -1) 
		{
			serror("fstat: %s", filename);
			(void) close(fd);
			continue;
		}

		/* These members, and seqno, vary over the loop. */
		status = set_timestamp(&prod.info.arrival);
		if(status != ENOERR)
			continue;

		if (key == (char *)NULL)
		{
			/* adding .* to the end of the filename for the ident */
			newid = (char *)calloc(1, (strlen(filename) + 3) * sizeof(char));
			if (newid == (char *)NULL)
			{
				exit(1);
			}
			strcat(newid, filename);
			strcat(newid, ".*");
			prod.info.ident = newid;
			free(newid);
		} else {
			prod.info.ident = key;
		}

		prod.info.sz = statb.st_size;
		prod.data = NULL;

#ifndef NO_MMAP
		prod.data = mmap(0, prod.info.sz,
			PROT_READ, MAP_PRIVATE, fd, 0);
		if(prod.data == NULL)
		{
			serror("mmap: %s", filename);
			(void) close(fd);
			continue;
		}

		if(mm_md5(md5ctxp, prod.data, prod.info.sz,
			prod.info.signature) != 0)
		{
			(void) munmap(prod.data, prod.info.sz);
			(void) close(fd);
			continue;
		}

		/*
		 * Do the deed
		 */
		status = pq_insert(pq, &prod);
#else /*!NO_MMAP*/

		if(fd_md5(md5ctxp, fd, statb.st_size, prod.info.signature) != 0)
		{
			(void) close(fd);
			continue;
		}
		if(lseek(fd, 0, SEEK_SET) == (off_t)-1)
		{
			serror("rewind: %s", filename);
			(void) close(fd);
			continue;
		}

		index = PQE_NONE;
		status = pqe_new(pq, &prod.info,
			&prod.data, &index);
		if(status != ENOERR)
			continue;
		if(read(fd, prod.data, prod.info.sz) !=
				prod.info.sz)
		{
			serror("read %s %u", filename, prod.info.sz);
			(void) pqe_discard(pq, index);
		}
		status = pqe_insert(pq, index);
#endif /*!NO_MMAP*/
		switch (status) {
		case ENOERR:
			/* no error */
			if(ulogIsVerbose())
				uinfo("%s",
	 s_prod_info(NULL, 0, &prod.info, ulogIsDebug())) ;
			break;
		case PQUEUE_DUP:
			uerror("Product already in queue: %s",
				s_prod_info(NULL, 0, &prod.info, 1));
			break;
		case ENOMEM:
			uerror("queue full?");
			break;	
		case EINTR:
#ifdef EDEADLOCK
		case EDEADLOCK:
#else
		case EDEADLK:
#endif
			/* TODO: retry ? */
		default:
			uerror("pq_insert: %s\n",
				status > 0 ? strerror(status) :
					"Internal error");
			break;
		}

#ifndef NO_MMAP
		(void) munmap(prod.data, prod.info.sz);
#endif /*!NO_MMAP*/
		(void) close(fd);
	}
	free_MD5_CTX(md5ctxp);	
	}

	exit(0);
}


