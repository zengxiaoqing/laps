#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "geoLib.h"

static FILE * bcd_file;
static int open_status = 0;

void init_bcd_file(char * bcd_file_name, int * status)
{
   char scratch[200];
   char * cp;
   int i;

   if (open_status) fclose(bcd_file);

   cp = getenv("GEO_DATA");
   if (cp!=NULL) {
       strcpy(scratch,cp);
       i = strlen(scratch);
       if (scratch[i-1]!='/') strcat(scratch,"/");
       strcat(scratch,bcd_file_name);
       bcd_file = fopen(scratch, "r");
   }/*endif*/
   if (bcd_file==NULL) bcd_file = fopen(bcd_file_name, "r");

   if (bcd_file==NULL) {
       open_status = 0;
       *status = 0;
   } else {
       *status = 1;
       open_status = 1;
   }/*endif*/
} /* end init_bcd_file */

void get_one_bcdrec(struct bcdRecord * bcd_rec, int *status)
{
   int n;

   if (open_status==0) {
       *status=0;
       return;
   }/*endif*/

   if (fread(&n,4,1,bcd_file)!=1 || n<2 || n>MAXPTS) {
       *status = 0;
       open_status = 0;
       return;
   }/*endif*/

   bcd_rec->npts = n;
   n = n*2+4;

   if (fread(&(bcd_rec->lat1),4,n,bcd_file)!=n) {
       *status = 0;
       open_status = 0;
       return;
   }/*endif*/

   *status = 1;

} /* end get_one_bcdrec */

