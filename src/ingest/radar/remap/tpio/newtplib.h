/*
 *	Newtplib.h...
 *
 *	Get the right tape interface.
 */

#ifndef	LINT
static char newtplib_h_id[] = "@(#)newtplib.h	2.1	2/15/95";
#endif	/* LINT */

#include "newtpdefs.h"

#ifdef	MT_NAMES
#include <sys/mtio.h>
#endif	/* MT_NAMES */

#ifdef	ST_NAMES
#include <sys/tape.h>

#define	MTBSF		STRSF
#define	MTBSR		STRSR
#define	MTFSF		STFSF
#define	MTFSR		STFSR
#define	MTIOCTOP	STIOCTOP
#define	MTREW		STREW
#define	MTWEOF		STWEOF
#define	MTOFFL		STOFFL
#define	mtop		stop
#define	mt_command	st_command
#define	mt_count	st_count
#define	mt_op		st_op

#endif	/* ST_NAMES */
