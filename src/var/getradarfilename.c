#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <dirent.h>

#ifdef FORTRANUNDERSCORE
#define getradarfilename getradarfilename_
#endif

#ifdef FORTRANDOUBLEUNDERSCORE
#define getradarfilename getradarfilename__
#endif

int getradarfilename(char *path,int *size,char *files,int *rtn)
{
   struct dirent  *dp;
   DIR            *dirp;

   *rtn=0;
   if ((dirp = opendir (path)) != NULL) {
      while ((dp = readdir (dirp)) != NULL) {
        if ( strncmp(dp->d_name ,".",1)==0 || strncmp(dp->d_name ,"..",2)==0 )
           continue;
        if ( strlen(files)+strlen(dp->d_name) > *size) {
           *rtn=-1;
           return(*rtn);
        }
        strcat(files,dp->d_name);
        strcat(files," ");
      }
      files[strlen(files)-1]='\0';
      *size=strlen(files);
      return(*rtn);
   }
}
