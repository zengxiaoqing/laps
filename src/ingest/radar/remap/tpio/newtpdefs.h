/*
 *	Newtpdefs.h...
 *
 *	Make some decisions and set some constants.
 */

#ifndef	LINT
static char newtpdefs_id[] = "@(#)newtpdefs.h	2.1	2/15/95";
#endif	/* LINT */

#define		MAX		20
#define		MAX_LEN		100
#define		MAX_ERRORS	20

#ifdef	SUN
#define	MT_NAMES
#endif

#ifdef	IBM
#define	ST_NAMES
#define	NO_UNDERSCORE
#endif

#ifdef	HP
#define	MT_NAMES
#endif

#ifdef	SGI
#define	MT_NAMES
#endif
