/*
 *	merge.h...
 *
 *	Info for the DBZ/VEL merger process.
 */

#ifdef	LINT
static char merge_id[] = "@(#)merge.h	2.2	2/14/95";
#endif	/* LINT */

#include <malloc.h>

#define	INCREMENT	16
#define	AZI_ACCEPT	100			/* one degree */
#define	MISSING_VALUE	0

#define	ALG_MISS	999.

/*
 *	Set minimum values and resolutions of the four fields.  All values
 *	are multiplied by 100.
 */

#define	DZ_MIN			-4500
#define	V_MIN			-6400
#define	SN_MIN			-3000
#define	SW_MIN			-3200

#define	DZ_RES			50
#define	V_RES			50
#define	SN_RES			50
#define	SW_RES			25

int rt_num_radials = 0;
int rt_current_alloc = 0;

int rt_using_dbz = 0;

struct dz_radial {
	short int first_gate;
	short int gate_spacing;
	short int number_of_gates;
	short int dummy;			/* for alignment */
	long int azimuth;
} *dz_radial;

int	dz_size = sizeof(struct dz_radial);

unsigned char **mptr;
