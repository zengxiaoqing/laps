/* new_fill.h */
#ifndef fill_laps_h
#define fill_laps_h

#define MAX_LAPS_NAME_LEN        10
#define LAPS_FN_LEN               9
#define MAX_WFO_VAR_NAME_LEN     20
#define MAX_WFO_LEVEL_NAME_LEN   12
#define MAX_WFO_LEVEL_LEN        5
#define MAX_WFO_LEVEL_VALUE      10
#define CHARS_PER_LEVEL          10
#define MAX_WFO_LEVELS           22
#define MAX_TIMESTRING_LEN       30
#define MAX_PARAMETER_LIST_LEN   100
#define MAX_FCST_HOURS 18
#define NUM_FCST_DAY 2
#define NO_MORE_FILES 5

#define TRUE    1
#define FALSE   0
#define SUCCESS TRUE
#define ERROR   FALSE

typedef struct  {
                char LAPS_dir_name    [MAX_LAPS_NAME_LEN];
                char LAPS_cdl_varname [MAX_LAPS_NAME_LEN];
                char LAPS_level_name  [MAX_LAPS_NAME_LEN];
                char WFO_cdl_varname  [MAX_WFO_VAR_NAME_LEN];
                char WFO_level_name   [MAX_WFO_LEVEL_NAME_LEN];
                char WFO_level        [MAX_WFO_LEVEL_LEN];
                char WFO_level_value  [MAX_WFO_LEVEL_VALUE];
                } PARAMETER_LIST_T;

typedef struct {
               char  WFO_level        [MAX_WFO_LEVEL_LEN];
               short WFO_level_value;
               } WFO_LEVELS_T;

int processLAPS(char *,time_t, time_t,int, PARAMETER_LIST_T *);
int fillWFOlevels(int, char *, WFO_LEVELS_T *);
int findWFOindex(int, WFO_LEVELS_T *, char *, char *, float *);
int extractLAPSdata(int *, char *, PARAMETER_LIST_T, char *,time_t, time_t, int *, int *, int *, float **, float **, short **);
int getDataSize(int, char *, int *, int *, int * );
int getLAPSdata(int, double, double, const char *,  int, int, int, float **, float **, short **);
int openOutputWFOfile(long *, time_t, long *, long *);
int transferLAPStoWFO(int, long, PARAMETER_LIST_T, time_t, time_t, int, int, int, float **, float **, short **);
int storeLAPSdata(long, int, double, double, char *, char *, char *, char *, int, int, int, float **, float **, short **);
int get_LAPS_cdf_filename (PARAMETER_LIST_T, char *,time_t, char *);
int get_WFO_cdf_filename (time_t, char *, char *);
int get_Laps_to_WFO_parameter_list (PARAMETER_LIST_T *, int *);
int copyTemplateToNewCDFfile(char *, char *);
void adj_neg_precip(int, int, float *);
void cvt_mm_to_meters(int, int, float *);
void cvt_in_to_meters(int, int, float *);
void filenameToUnixtime(time_t *, char *);
void jdayToMoDay(int, int, int *, int *);
long dayInfo2unixtime(int, int, int, int, int, int);
int incrementFcstHour(char *, char *, long *);
int incrementFcstPeriod(time_t, char *);
			
#endif  /* Don't put anything after this line! */

