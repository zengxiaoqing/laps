/*
 *	defs.h...
 *
 *	Definitions...
 */

#ifdef	LINT
static char defs_id[] = "@(#)defs.h	2.2	2/14/95";
#endif	/* LINT */

#define	DBZ			0
#define	VEL			1
#define	SPW			2
#define	SNR			3
#define	CMP			4
#define	DBM			5	/* not currently used */
#define	PCP1			6	/* created in RT_PRECIP */
#define	PCP3			7	/* created in RT_PRECIP */
#define	PCPT			8	/* created in RT_PRECIP */
#define	RVEL			9	/* not currently used */
#define	ZDR			10	/* differential reflectivity */
#define	PDP			11	/* differential phase */
#define	RHV			12	/* differential correlation coef */
#define	KDP			13	/* ??? not currently used */

#define	NON_EXISTENT		-888	/* will phase this out */
#define	MISSING			-999

#define	MISSING_DATA		-999.0
#define	RANGE_FOLDED_DATA	-888.0
