typedef struct CDF_GRID_INFO {
	int varid;
	int level_coord;
	int fctime_coord;
	int	qptr;
	int x_dim;
	int y_dim;
}cdf_grid_info;

typedef struct GRID_ID_INFO {
	char model[10];
	char field[10];
	char units[10];
	char level[10];
	int  fctime;
	int  i4time;
}grid_id_info;

#ifdef __alpha
static int nbytes[]={0,1,1,2,8,4,8};	/* #bytes in each NC TYPE */
#else
static int nbytes[]={0,1,1,2,4,4,8};	/* #bytes in each NC TYPE */
#endif

static int N_headers = 0;
