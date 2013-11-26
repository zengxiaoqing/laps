#include <dirent.h>
#include <stdio.h>
#include <sys/shm.h>
#include <fcntl.h>
#include <sys/stat.h>

#ifdef FORTRANUNDERSCORE
#define getfilenames getfilenames_
#endif

#ifdef FORTRANDOUBLEUNDERSCORE
#define getfilenames getfilenames__
#endif

int getfilenames(char *path,int *size,char *files,int *rtn)
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
/* added: return here in response to compiler warning
 * getiofile.c(36): warning #1011: missing return statement 
 * at end of non-void function "getfilenames_". Hongli Jiang 10/21/2013 */
     return(*rtn);
}


#ifdef FORTRANUNDERSCORE
#define rioreadfile rioreadfile_
#endif

#ifdef FORTRANDOUBLEUNDERSCORE
#define rioreadfile rioreadfile__
#endif

int rioreadfile(char *fname,char *buf,int *size,int *rtn_code)
{
      struct stat st_file;
      int    fd_file;

     /** check if file exist, if not, return error */
     if (stat(fname, &st_file) != 0){
        printf("file %s not found\n",fname);
        *rtn_code=-1;
        return(*rtn_code);
     }

     /** check if buffer size smaller than file size, if yes,
         return error */
     if ( st_file.st_size > *size ) {
        printf("buffer size too small,real size is %ld",st_file.st_size);
        *rtn_code=-2;
        return(*rtn_code);
     }

     /** open file */
     if ((fd_file = open(fname,O_RDONLY)) < 0){
        *rtn_code=-1;
        return(*rtn_code);
     }

     /** read data from file and put into buffer */
     *size=read(fd_file,buf,st_file.st_size);
     *rtn_code=0;
     close(fd_file);
     return(*rtn_code);
}
