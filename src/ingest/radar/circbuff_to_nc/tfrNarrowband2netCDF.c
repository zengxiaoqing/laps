
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <dirent.h>
#include <unistd.h>
#include <netcdf.h>
#include "RadarAccessor.H"

#define TRUE 0
#define FALSE 1
#define ERROR -1
#define FILE_EXISTS -2
#define SUCCESS 0
#define MAX_FN_LEN 200
#define MAX_TIMESTRING_LEN 20
#define SYSCMD "ncgen -o %s %s"
#define FILENAME_LEN 14
#define RAD_SITE_NAME_LEN 132
#define RAD_MNEMONIC_LEN 5
#define RADAR_INFO_SITE "/awips/fxa/data/localization/nationalData/radarInfoMaster.txt"

int scaleN2W(RadialReport *);
int writeNarrowband(RadialReport *, char *, float, float, char *, time_t);
void free_malloc(char *, char *, char *, RadialReport *);
void downcase_c(char *,char *);
int filename2Unixtime(time_t *, char *);
long dayInfo2unixtime(int, int, int, int, int, int);
int getOutputFilename (char *, time_t, char *, int *);
int openOutputFile (char *);
int getSiteInfo(char *, float *, float *, float *);
bool decodeRadial(TextString, float, RadialReport *, bool);

// NOTE: The environment variable $LAPS_DATA_ROOT needs to be defined before
//       running this program

//pass in homeRadar startTimeWindow endTimeWindow inputDirName outputDir
// homeRadar should be 4 characters, e.g. "kftg"
// startTimeWindow and endTimeWindow are time_t and are in unixtime
// inputDirName and outputDir are character strings

int main (int argc, char *argv[])
    {
	DIR *inputDir;
	char homeRadar[RAD_MNEMONIC_LEN], fileName[FILENAME_LEN]; 
	char *inputDirName, *outputDir, *fullFileName;
	bool ring, b_status;
	int status, productCode, level, count_proc; 
	float elevAngle, resolution;
	time_t startTimeWindow, endTimeWindow, fileTime;
	RadialReport *report;
	struct dirent dir_ent, *p_dirent;

        if (argc != 6) {  // Can't run...don't have all the parameters you need
	  printf("xferNarrowband program needs to be called with 5 parameters.\n");
	  printf("Look at calling script for requirements.\n");
	  exit(ERROR);
	}

// capture the arguments from the command line
	strcpy(homeRadar, argv[1]);
	downcase_c(homeRadar, homeRadar);
	startTimeWindow = (time_t)atoi(argv[2]);
	endTimeWindow = (time_t)atoi(argv[3]);
	inputDirName = (char *)malloc((strlen(argv[4]) + 5) * sizeof(char));
	strcpy(inputDirName, argv[4]);
	outputDir = (char *)malloc((strlen(argv[5]) + 5) * sizeof(char));
	strcpy(outputDir, argv[5]);

// setup fullFileName and report for use
	fullFileName = (char *)malloc(MAX_FN_LEN * sizeof(char));
	report = (RadialReport *)malloc(sizeof(RadialReport));

// no "missing data" ring around outside of data
	ring = false; 

// for right now retrieving only product code 20, resol 2, elev 0.5, level 16
	productCode = 20;
	resolution = 2.0;  
	level = 16;
	elevAngle = 0.5;

// open inputDir, find files within TimeWindow, and process them
	inputDir = opendir(inputDirName);
        if (inputDir == NULL) {
          fprintf(stdout, "Unable to open directory %s.\n", inputDirName);
	  free_malloc(inputDirName, outputDir, fullFileName, report);
 	  exit(ERROR);
        }

// expects filename of the format YYYYMMDD_HHMM
	p_dirent = &dir_ent;
        p_dirent = readdir(inputDir);
        fileName[13] = '\0';
	count_proc = 0;
        while (p_dirent != NULL) {
	  strcpy(fileName,p_dirent->d_name);
          if (strncmp(fileName,".",1) == 0) {
// don't process . and ..
          } 
          else {
	    status = filename2Unixtime(&fileTime, fileName);
	    if (status == ERROR) {
	      printf("ERROR converting filename %s to unixtime.\n", fileName);
	    }
	    else {
	      if ((fileTime >= startTimeWindow) && (fileTime <= endTimeWindow)) {
	        count_proc += 1;  // keeps track of number of files in time window

// make fullFileName which has directory and filename
		strcpy(fullFileName, inputDirName);
		strcat(fullFileName, "/");
		strcat(fullFileName, fileName);

// read radial data from file and return in report
	        b_status = decodeRadial((TextString)fullFileName, 
                                        resolution, report, ring);

// scale data to wideband counts
	   	status = scaleN2W(report);

// write out data in report
	        status = writeNarrowband(report, homeRadar, elevAngle,
                                         resolution, outputDir, fileTime);
		if (status != SUCCESS) {
		  if (status == FILE_EXISTS) printf("Return status from writeNarrowband: %d FILE EXISTS\n",status);
		  if (status == ERROR) printf("Return status from writeNarrowband: %d ERROR\n",status);
		}
		delete report->radialAngles;
		delete report->radialData;
	      }
	    }
          }
          p_dirent = readdir(inputDir);

	}

	closedir(inputDir);
	if (count_proc == 0) printf("No radar data files available within time window.\n");
	free_malloc(inputDirName, outputDir, fullFileName, report);
	exit(SUCCESS);
    }
/**************************************************************************************/
int scaleN2W(RadialReport *report)
    {
	int i, gatesXrad;
	long VCP;
	unsigned char *dptr;
	unsigned char clear[16], storm[16];

	clear[0] = 0;
	clear[1] = 10; 
	clear[2] = 18; 
	clear[3] = 26; 
	clear[4] = 34; 
	clear[5] = 42; 
	clear[6] = 50; 
	clear[7] = 58; 
	clear[8] = 66; 
	clear[9] = 74; 
	clear[10] = 82; 
	clear[11] = 90; 
	clear[12] = 98; 
	clear[13] = 106; 
	clear[14] = 114; 
	clear[15] = 122; 

	storm[0] = 0;
	storm[1] = 76;
	storm[2] = 86;
	storm[3] = 96;
	storm[4] = 106;
	storm[5] = 116;
	storm[6] = 126;
	storm[7] = 136;
	storm[8] = 146;
	storm[9] = 156;
	storm[10] = 166;
	storm[11] = 176;
	storm[12] = 186;
	storm[13] = 196;
	storm[14] = 206;
	storm[15] = 216;

	VCP = report->VCP;
	gatesXrad = report->gatesPerRadial * report->numRadials;

	dptr = (unsigned char *)report->radialData;

	if (report->mode == ClearAir) {
	   printf("VCP = %d...mode = %d...ClearAir\n",report->VCP, report->mode);
	}
	if (report->mode == SevereWx) {
	  printf("VCP = %d...mode = %d SevereWx\n",report->VCP, report->mode);
	}

	if (report->mode == Maintenance) {
	  printf("Radar in Maintenance mode...no scaling.\n");
	}
	else if ((VCP == 31 || VCP == 32) && (report->mode == ClearAir)) {
	  printf("Using Clear Air scaling.\n");
	  for (i = 0; i < gatesXrad; i++) {
	    *dptr = clear[((*dptr >> 4) & 0x0f)];
	    dptr++;
	  }
	}
	else {
	  printf("Using Storm scaling.\n");
	  for (i = 0; i < gatesXrad; i++) {
	    *dptr = storm[((*dptr >> 4) & 0x0f)];
	    dptr++;
	  }
	}

	return(SUCCESS);

    }
/**************************************************************************************/
int writeNarrowband(RadialReport *report, char *homeRadar, float elevAngle,
                    float resolution, char *outputDir, time_t fileTime)
    {

	short shortVal, numGatesV;
	int cdfid, varid, dimid, elevNum, status, i, doubleLoc;
	int radarNameLen_file, radarNameLen_var; 
	float gateSizeZ, gateSizeV, lat, lon, siteElev;
        float firstGateRangeZ, firstGateRangeV;
	double doubleVal;
	size_t dim_len, mindex[1], start_1[1], count_1[1], start_2[2], count_2[2];
	unsigned char *radialData;
	char fileName[MAX_FN_LEN], cdlFile[MAX_FN_LEN]; 

    //  setup ancillary data for output to file
        numGatesV = 920;
        firstGateRangeZ = 0;
        firstGateRangeV = -0.375;
        gateSizeV = .25;  // in km

    //  Make output filename: format YYJJJHHMM_elevXX, where XX is elevation number
	status = getOutputFilename(outputDir, fileTime, fileName, &elevNum);

    //  Open netCDF output file
	cdfid = openOutputFile(fileName);
	if (cdfid == FILE_EXISTS) {
	  return(FILE_EXISTS);
	}
	else if (cdfid == ERROR) {
	  return(ERROR);
	}

    //  Setup for netCDF write
	mindex[0] = 0;

	start_1[0] = 0;

	start_2[0] = 0;
	start_2[1] = 0;

    //  Verify Z_bin dimension >= report->gatesPerRadial
	status = nc_inq_dimid(cdfid,"Z_bin", &dimid);
	if (status != NC_NOERR) {
	  printf("Cannot determine number of gates output file can hold.\n");
          status = nc_close(cdfid);
          return(ERROR);
        }

        status = nc_inq_dimlen(cdfid, dimid, &dim_len);
        if (status != NC_NOERR) {
	  printf("Cannot determine number of gates output file can hold.\n");
          status = nc_close(cdfid);
          return(ERROR);
        }

	if ((int)dim_len < report->gatesPerRadial) {
	  printf("Output file %s can only hold %d gates of radar data.\n", fileName);
	  printf("Radar report has %d gates of radar data.\n", report->gatesPerRadial);
          status = nc_close(cdfid);
	  return(ERROR);
	}

    // if only 230 gates, fill 460 gates in file by doubling data from 230 gates
	if (report->gatesPerRadial == 230) {
	  radialData = (unsigned char *)malloc(460 * report->numRadials);
	  doubleLoc = 0;
	  for (i = 0; i < (report->gatesPerRadial*report->numRadials); i++) {
	    radialData[doubleLoc] = (unsigned char)report->radialData[i];
	    radialData[doubleLoc+1] = (unsigned char)report->radialData[i];
	    doubleLoc += 2;
	  }
	}
	else {
	  radialData = (unsigned char *)malloc(report->gatesPerRadial * report->numRadials);
	  for (i = 0; i < (report->gatesPerRadial*report->numRadials); i++) {
	    radialData[i] = (unsigned char)report->radialData[i];
	  }
	}

    //  Write elevationNumber (calculated in getOutputFilename)
        status = nc_inq_varid(cdfid, "elevationNumber", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'elevationNumber' not found in output file.\n");
          status = nc_close(cdfid);
	  free(radialData);
	  return(ERROR);
        }

	shortVal = (short)elevNum;
        status = nc_put_var1_short(cdfid, varid, mindex, &shortVal);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'elevNum' to output file %s.\n", fileName);
          status = nc_close(cdfid);
	  free(radialData);
          return(ERROR);
        }

    //  Write elevationAngle from elevAngle
        status = nc_inq_varid(cdfid, "elevationAngle", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'elevationAngle' not found in output file.\n");
          status = nc_close(cdfid);
	  free(radialData);
	  return(ERROR);
        }

        status = nc_put_var1_float(cdfid, varid, mindex, &elevAngle);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'elevAngle' to output file %s.\n", fileName);
          status = nc_close(cdfid);
	  free(radialData);
          return(ERROR);
        }

    //  Write numRadials from report->numRadials
        status = nc_inq_varid(cdfid, "numRadials", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'numRadials' not found in output file.\n");
          status = nc_close(cdfid);
	  free(radialData);
	  return(ERROR);
        }

	shortVal = (short)report->numRadials;
        status = nc_put_var1_short(cdfid, varid, mindex, &shortVal);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'report->numRadials' to output file %s.\n", fileName);
          status = nc_close(cdfid);
	  free(radialData);
          return(ERROR);
        }

    //  Write radialAzim from report->radialAngles
        status = nc_inq_varid(cdfid, "radialAzim", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'radialAzim' not found in output file.\n");
          status = nc_close(cdfid);
	  free(radialData);
	  return(ERROR);
        }

	count_1[0] = report->numRadials;
        status = nc_put_vara_float(cdfid, varid, start_1, count_1, report->radialAngles);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'report->radialAngles' to output file %s.\n", fileName);
          status = nc_close(cdfid);
	  free(radialData);
          return(ERROR);
        }

    //  Write Z (radial, num_gates) from radialData
        status = nc_inq_varid(cdfid, "Z", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'Z' not found in output file.\n");
          status = nc_close(cdfid);
	  free(radialData);
	  return(ERROR);
        }

	count_2[0] = report->numRadials;
	if (report->gatesPerRadial == 230) {
	  count_2[1] = 460;
	}
	else {
	  count_2[1] = report->gatesPerRadial;
	}

        status = nc_put_vara_uchar(cdfid, varid, start_2, count_2, radialData);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'radialData' to output file %s.\n", fileName);
          status = nc_close(cdfid);
	  free(radialData);
          return(ERROR);
        }
	free(radialData);

    //  get site info from LAPS external file (using homeRadar to get lat, 
    //      lon, siteElev from "/awips/fxa/data/localization/nationalData/radarInfoMaster.txt")
	status = getSiteInfo(homeRadar, &lat, &lon, &siteElev);
	if (status == ERROR) {
	  printf("ERROR getting radar site info from RADAR_INFO_SITE\n");
          status = nc_close(cdfid);
          return(ERROR);
	}

    //  Write radarName
	status = nc_inq_dimid(cdfid,"radarNameLen", &dimid);
	if (status != NC_NOERR) {
	  printf("Cannot determine dimension 'radarNameLen' in output file.\n");
          status = nc_close(cdfid);
          return(ERROR);
        }

        status = nc_inq_dimlen(cdfid, dimid, &dim_len);
        if (status != NC_NOERR) {
	  printf("Cannot determine dimension 'radarNameLen' in output file.\n");
          status = nc_close(cdfid);
          return(ERROR);
        }
	radarNameLen_file = (int)dim_len;
	radarNameLen_var = strlen(homeRadar);

	if (radarNameLen_var > radarNameLen_file) {
	  count_1[0] = radarNameLen_file - 1;
	}
	else {
	  count_1[0] = radarNameLen_var;
	}

        status = nc_inq_varid(cdfid, "radarName", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'radarName' not found in output file.\n");
          status = nc_close(cdfid);
	  return(ERROR);
        }

        status = nc_put_vara_text(cdfid, varid, start_1, count_1, homeRadar);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'homeRadar' to output file %s.\n", fileName);
          status = nc_close(cdfid);
          return(ERROR);
        }

    //  Write siteLat
        status = nc_inq_varid(cdfid, "siteLat", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'siteLat' not found in output file.\n");
          status = nc_close(cdfid);
	  return(ERROR);
        }

        status = nc_put_var1_float(cdfid, varid, mindex, &lat);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'lat' to output file %s.\n", fileName);
          status = nc_close(cdfid);
          return(ERROR);
        }

    //  Write siteLon
        status = nc_inq_varid(cdfid, "siteLon", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'siteLon' not found in output file.\n");
          status = nc_close(cdfid);
	  return(ERROR);
        }

        status = nc_put_var1_float(cdfid, varid, mindex, &lon);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'lon' to output file %s.\n", fileName);
          status = nc_close(cdfid);
          return(ERROR);
        }

    //  Write siteAlt
        status = nc_inq_varid(cdfid, "siteAlt", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'siteAlt' not found in output file.\n");
          status = nc_close(cdfid);
	  return(ERROR);
        }

        status = nc_put_var1_float(cdfid, varid, mindex, &siteElev);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'siteElev' to output file %s.\n", fileName);
          status = nc_close(cdfid);
          return(ERROR);
        }

    //  Write VCP (volume coverage pattern) 
        status = nc_inq_varid(cdfid, "VCP", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'VCP' not found in output file.\n");
          status = nc_close(cdfid);
	  return(ERROR);
        }

	shortVal = (short)report->VCP;
        status = nc_put_var1_short(cdfid, varid, mindex, &shortVal);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'report->VCP' to output file %s.\n", fileName);
          status = nc_close(cdfid);
          return(ERROR);
        }

    //  Write esEndTime from report->obsTime
        status = nc_inq_varid(cdfid, "esEndTime", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'esEndTime' not found in output file.\n");
          status = nc_close(cdfid);
	  return(ERROR);
        }

	doubleVal = (double)fileTime;
        status = nc_put_var1_double(cdfid, varid, mindex, &doubleVal);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'fileTime' to output file %s.\n", fileName);
          status = nc_close(cdfid);
          return(ERROR);
        }

    //  Write firstGateRangeZ 
        status = nc_inq_varid(cdfid, "firstGateRangeZ", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'firstGateRangeZ' not found in output file.\n");
          status = nc_close(cdfid);
	  return(ERROR);
        }

        status = nc_put_var1_float(cdfid, varid, mindex, &firstGateRangeZ);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'firstGateRangeZ' to output file %s.\n", fileName);
          status = nc_close(cdfid);
          return(ERROR);
        }

    //  Write firstGateRangeV 
        status = nc_inq_varid(cdfid, "firstGateRangeV", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'firstGateRangeV' not found in output file.\n");
          status = nc_close(cdfid);
	  return(ERROR);
        }

        status = nc_put_var1_float(cdfid, varid, mindex, &firstGateRangeV);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'firstGateRangeV' to output file %s.\n", fileName);
          status = nc_close(cdfid);
          return(ERROR);
        }

    //  Write gateSizeZ from resolution
	if ((resolution == 1.0) || (report->gatesPerRadial == 230)) {
	  gateSizeZ = 1.0;  // .54 nm = 1.0 km
	}
	else if (resolution == 2.0) {
	  gateSizeZ = 2.0;  // 1.08 nm = 2.0 km
	}
	else {
	  printf("Resolution of %f is unknown...should be 1.0 or 2.0.\n", resolution);
          status = nc_close(cdfid);
	  return(ERROR);
	}

        status = nc_inq_varid(cdfid, "gateSizeZ", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'gateSizeZ' not found in output file.\n");
          status = nc_close(cdfid);
	  return(ERROR);
        }

        status = nc_put_var1_float(cdfid, varid, mindex, &gateSizeZ);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'gateSizeZ' to output file %s.\n", fileName);
          status = nc_close(cdfid);
          return(ERROR);
        }

    //  Write gateSizeV 
        status = nc_inq_varid(cdfid, "gateSizeV", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'gateSizeV' not found in output file.\n");
          status = nc_close(cdfid);
	  return(ERROR);
        }

        status = nc_put_var1_float(cdfid, varid, mindex, &gateSizeV);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'gateSizeV' to output file %s.\n", fileName);
          status = nc_close(cdfid);
          return(ERROR);
        }

    //  Write numGatesZ from report->gatesPerRadial
	if (report->gatesPerRadial == 230) {
	  shortVal = 460;
	}
	else {
	  shortVal = (short)report->gatesPerRadial;
	}

        status = nc_inq_varid(cdfid, "numGatesZ", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'numGatesZ' not found in output file.\n");
          status = nc_close(cdfid);
	  return(ERROR);
        }

        status = nc_put_var1_short(cdfid, varid, mindex, &shortVal);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'report->gatesPerRadial' to output file %s.\n", fileName);
          status = nc_close(cdfid);
          return(ERROR);
        }

    //  Write numGatesV 
        status = nc_inq_varid(cdfid, "numGatesV", &varid);
        if (status != NC_NOERR) {
          printf("Variable 'numGatesV' not found in output file.\n");
          status = nc_close(cdfid);
	  return(ERROR);
        }

        status = nc_put_var1_short(cdfid, varid, mindex, &numGatesV);
        if (status != NC_NOERR) {
	  printf("ERROR writing 'numGatesV' to output file %s.\n", fileName);
          status = nc_close(cdfid);
          return(ERROR);
        }

    //  Close netCDF output file
        status = nc_close(cdfid);
        if (status != NC_NOERR) {
          printf("ERROR closing netCDF file %s.\n", fileName);
          return(ERROR);
        }
	else {
	  printf("Radar file %s successfully written.\n", fileName);
	  return(SUCCESS);
	}

    }
/**************************************************************************************/
void free_malloc(char *inputDirName, char *outputDir, char *fullFileName,
		 RadialReport *report)
{
	free(inputDirName);
	free(outputDir);
	free(fullFileName);
	free(report);
	return;
}
/**************************************************************************************/
void downcase_c(char *instr,char *outstr)
{
        char *inptr, *tempout;
        short i,len_in;
        char outchar[2];

        outchar[1] = '\0';
        len_in = strlen(instr);
        tempout = (char *) malloc(len_in + 1);
        tempout[0] = '\0';
        tempout[1] = '\0';
        inptr = instr;

        for (i = 0; i < len_in; i++) {
          outchar[0] = tolower(*inptr);
          if (i == 0)
            strncpy(tempout,outchar,1);
          else
            strncat(tempout,outchar,1);
          inptr++;
        }

        strcpy(outstr,tempout);
        free(tempout);
}
/**************************************************************************************/
int filename2Unixtime(time_t *refTime, char *filename)
{
        int hour_r, min_r, yr, mo, day, sec, fn_len;
        char c_hour_r[3], c_min_r[3];
        char c_yr[3], *p, c_mo[3], c_day[3];
        char timeString[30];

        fn_len = strlen(filename);
        p = filename + 8;
        if (fn_len == 13 && strncmp(p,"_",1) == 0) { /* wfo filename format */
          p = filename + 2;
          strncpy(c_yr,p,2);
          c_yr[2] = '\0';
          yr = atoi(c_yr);

          p = filename + 4;
          strncpy(c_mo,p,2);
          c_mo[2] = '\0';
          mo = atoi(c_mo);

          p = filename + 6;
          strncpy(c_day,p,2);
          c_day[2] = '\0';
          day = atoi(c_day);

          p = filename + 9;
          strncpy(c_hour_r,p,2);
          c_hour_r[2] = '\0';
          hour_r = atoi(c_hour_r);

          p = filename + 11;
          strncpy(c_min_r,p,2);
          c_min_r[2] = '\0';
          min_r = atoi(c_min_r);

        }
	else {
	  printf("Filename: %s is not of format YYYYMMDD_HHMM.\n", filename);
	  return(ERROR);
        }

        sec = 0;
        *refTime = (time_t)dayInfo2unixtime(yr, mo, day, hour_r, min_r, sec);

        return(SUCCESS);

}
/**************************************************************************************/
/* Y2K note:  This program will only work until 12/31/2069.  After that
              it will overlap with the 1/1/1970 start date of unixtime  */

long dayInfo2unixtime(int yr, int mo, int day, int hour, int min, int sec)
{

        int dayInMo[12], n_leap_1970_to_1999;
        long nsecmo[12];

        long nsecyr = 31536000, nsecda = 86400, nsechr = 3600, nsecmn = 60;
        long sum, nyr, nyrs, nleap, prevLeapYear, ibase, i;

        dayInMo[0] = 31;
        dayInMo[1] = 28;
        dayInMo[2] = 31;
        dayInMo[3] = 30;
        dayInMo[4] = 31;
        dayInMo[5] = 30;
        dayInMo[6] = 31;
        dayInMo[7] = 31;
        dayInMo[8] = 30;
        dayInMo[9] = 31;
        dayInMo[10] = 30;
        dayInMo[11] = 31;

        nsecmo[0] = 2678400;
        nsecmo[1] = 2419200;
        nsecmo[2] = 2678400;
        nsecmo[3] = 2592000;
        nsecmo[4] = 2678400;
        nsecmo[5] = 2592000;
        nsecmo[6] = 2678400;
        nsecmo[7] = 2678400;
        nsecmo[8] = 2592000;
        nsecmo[9] = 2678400;
        nsecmo[10] = 2592000;
        nsecmo[11] = 2678400;

        ibase = 70;
        prevLeapYear = 68;
        n_leap_1970_to_1999 = 7;

/* sum the number of years */
        nyr = yr;
        if (nyr > 1900) nyr = nyr - 1900;
        nyrs = nyr - ibase;

        if (nyr >= 70 && nyr <= 99) {
          nleap = ((nyr - prevLeapYear) / 4);
        }
        else if (nyr >= 0 && nyr < 70){
          nyrs = nyr + 30;
          prevLeapYear = 96;
          nleap = ((nyr + 100 - prevLeapYear) / 4) + n_leap_1970_to_1999;
        }
        else { /* date is before 1/1/1970 or after 12/31/2069 */
          fprintf(stdout,"Date is before 1/1/1970 or after 12/31/2069.\n");
          fprintf(stdout,"Function filename2Unixtime will not work.\n");
          return ERROR;
        }
        sum = nyrs * nsecyr + (nleap * nsecda);

/* sum the number of months */
        if ((mo < 1) || (mo > 12)) return ERROR;
        if (mo != 1) {
          for (i = 0; i < (mo - 1); i++) {
            sum = sum + nsecmo[i];
          }
        }

/* correct for Jan of Feb of a leap year */

        if ((nyr % 4) == 0) {
          if (mo <= 2) sum = sum - nsecda;
        }
/* sum the number of days */
        if ((day < 1) || (day > dayInMo[mo - 1])) return(ERROR);
        sum = sum + ((day - 1) * nsecda);

/* sum the number of hours */
        if ((hour < 0) || (hour > 23)) return(ERROR);
        sum = sum + (hour * nsechr);

/* sum the number of mins */
        if ((min < 0) || (min > 59)) return(ERROR);
        sum = sum + (min * nsecmn);

/* sum the number of secs */
        if ((sec < 0) || (sec > 59)) return(ERROR);
        sum = sum + sec;

        return sum;
}
/**************************************************************************************/
int getOutputFilename (char *outputDir, time_t fileTime, char *fileName, int *i_elevNum)
{

	char timeString[MAX_TIMESTRING_LEN], elevNum[3];
	struct tm *time = gmtime (&fileTime);

    // Make output filename: format YYJJJHHMM_elevXX, where XX is elevation number
	strftime(timeString, MAX_TIMESTRING_LEN, "%y%j%H%M", time);
	strcpy(fileName, outputDir);
	strcat(fileName, "/");
	strcat(fileName, timeString);
	strcat(fileName, "_elev");

    // Look in output directory to get elevationNumber 
	*i_elevNum = 1;
	strcpy(elevNum, "01");

    // Add _elev## to filename 
	strcat(fileName, elevNum);

	return(SUCCESS);
}
/**************************************************************************************/
int openOutputFile (char *outputFile)
{
	int outLen, cdlLen, cdfid, status, nc_noerr;
	char *syscmd, cdlFile[MAX_FN_LEN];	

    //  See if file is already there
	if( access(outputFile, F_OK) != 0 ) { // file does not exist 

    //  Ncgen netCDF output file (from $LAPS_DATA_ROOT/cdl/narrowband.cdl)
	  strcpy(cdlFile, getenv ("LAPS_DATA_ROOT"));
	  strcat(cdlFile, "/cdl/narrowband.cdl");
	  cdlLen = strlen(cdlFile);
	  outLen = strlen(outputFile);
	  syscmd = (char *)malloc((strlen(SYSCMD) + outLen + cdlLen + 10) * sizeof(char));
          sprintf(syscmd,SYSCMD, outputFile, cdlFile);

    //  create file, then open it 
          system(syscmd);
	  free(syscmd);

	  nc_noerr = NC_NOERR;
          status = nc_open(outputFile,NC_WRITE, &cdfid);
          if (status == NC_NOERR) {
            return(cdfid);
	  }
	  else {
	    printf("Could not open output file %s.\n", outputFile);
	    return(ERROR);
          }
        }
	else {  // file is already there...don't process again
	  printf("Radar file %s already processed.\n", outputFile);
	  return(FILE_EXISTS);
	}
}
/**************************************************************************************/
int getSiteInfo(char *radarMnemonic, float *lat, float *lon, 
                float *siteElev)
{

	FILE *fp;
	int found, radLen;
	char radarFile[MAX_FN_LEN], radarName[RAD_MNEMONIC_LEN], paramString[45];
	char junk1[4], junk2[5], junk3[5], s_lat[9], s_lon[11], s_elev[9];

	radLen = RAD_MNEMONIC_LEN - 1;
	strcpy(radarFile, RADAR_INFO_SITE);

	if( access(radarFile, F_OK) != 0 ) { // file does not exist 
          printf ("Radar site info file %s does not exist.\n", radarFile);
          return(ERROR);
        }

	if ((fp=fopen (radarFile, "r")) == (FILE *) 0) {
          printf ("Unable to open %s\n", radarFile);
          return(ERROR);
        }

	found = FALSE;
	while (!(feof (fp)) && found == FALSE) {
          fgets (paramString, 45, fp);
          if (strlen (paramString) > 0  && paramString[0] != '#') {
            sscanf (paramString, "%s %s %s %s %s %s %s", 
                                 radarName,
				 s_lat,
				 s_lon,	
				 junk1,
				 junk2,
				 junk3,
				 s_elev);
	    *lat = atof(s_lat);
	    *lon = atof(s_lon);
	    *siteElev = atof(s_elev);
	    if (strncmp(radarName, radarMnemonic,radLen) == 0) found = TRUE;
	  }
	}
	
	fclose(fp);

	if (found == TRUE) {
	  return(SUCCESS);
	}
	else {
	  printf("ERROR finding radar %s in %s.\n", *radarMnemonic, radarFile);
	  return(ERROR);
	}
}

