/*
 *	Newtpcfor.h...
 *
 *	Remap fortran interface names, if necessary.
 */

#ifndef	LINT
static char newtpcfor_h_id[] = "@(#)newtpcfor.h	2.1	2/15/95";
#endif	/* LINT */

#include "newtpdefs.h"

#ifdef	NO_UNDERSCORE

#define	tpopen_		tpopen
#define	tpread_		tpread
#define	tpwrite_	tpwrite
#define	tpclose_	tpclose
#define	tprew_		tprew
#define	tpbsf_		tpbsf
#define	tpbsr_		tpbsr
#define	tpfsf_		tpfsf
#define	tpfsr_		tpfsr
#define	tpweof_		tpweof
#define	tpweoff_	tpweoff
#define	tpstatus_	tpstatus
#define	tpseek_		tpseek

#endif	/* NO_UNDERSCORE */
