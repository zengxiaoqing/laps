#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef FORTRANUNDERSCORE
#define read_gms_image read_gms_image_
#endif
#ifdef FORTRANDOUBLEUNDERSCORE
#define read_gms_image read_gms_image__
#endif
#ifdef FORTRANCAPS
#define read_gms_image READ_GMS_IMAGE
#endif

#define Pi 3.141592654
#define SIZE 100
#define Width 1 
#define FILE_OUT_PIXEL "pixel"
#define DIMEN 1024
#define DIMEN_CAL 256

float  mercatorProjLat(float latitude);
float  inverseMercatorProjLat(float pt);

#ifdef WORDS_BIGENDIAN
#else
void swapShortInt(short int *n);
void swapInt(int *n);
void swapFloat(float *n);
#endif

void read_gms_image (char *filename,  float image_data[DIMEN][DIMEN], float calibration_table[DIMEN_CAL], int* istat ) {

	char *satfile;

   	FILE *fp,*ft;
	short int channel, no_grid_pt;
	int width_image_pixel,height_image_line;
	int offset_cal, offset_imageData, offset_gridData;

	float grid_lat[SIZE], grid_lon[SIZE]; 
        short int pixel_number_of_centre[SIZE],  line_number_of_centre[SIZE];
	float lat, lon, inverselat;
	int i=0,j, offset;
	unsigned char pixel_value;
	float	cal_value;
	int missing_count=0;

	satfile = filename;

	fp=fopen(satfile, "rb");
	if (fp == NULL) *istat = 0;   /* fail opening file */
/*
	ft=fopen(FILE_OUT_PIXEL, "w");
*/

        fseek(fp,12,SEEK_SET);
        fread(&channel,2,1,fp);
#ifdef WORDS_BIGENDIAN
#else
        swapShortInt(&channel);
#endif

        fseek(fp,36,SEEK_SET);
        fread(&width_image_pixel,4,1,fp);
#ifdef WORDS_BIGENDIAN
#else
        swapInt(&width_image_pixel);
#endif

        fseek(fp,40,SEEK_SET);
        fread(&height_image_line,4,1,fp);
#ifdef WORDS_BIGENDIAN
#else
        swapInt(&height_image_line);
#endif

        fseek(fp,72,SEEK_SET);
        fread(&offset_cal,4,1,fp);
#ifdef WORDS_BIGENDIAN
#else
        swapInt(&offset_cal);
#endif

        fseek(fp,76,SEEK_SET);
        fread(&offset_imageData,4,1,fp);
#ifdef WORDS_BIGENDIAN
#else
        swapInt(&offset_imageData);
#endif

        fseek(fp,80,SEEK_SET);
        fread(&no_grid_pt,2,1,fp);
#ifdef WORDS_BIGENDIAN
#else
        swapShortInt(&no_grid_pt);
#endif

        fseek(fp,84,SEEK_SET);
        fread(&offset_gridData,4,1,fp);
#ifdef WORDS_BIGENDIAN
#else
        swapInt(&offset_gridData);
#endif

/*
          fprintf(ft,"%s\n",satfile);
*/
    
	i=0;
        for (offset=offset_gridData; offset<offset_gridData+no_grid_pt*28; offset+=28){


		fseek(fp,offset,SEEK_SET);
		fread(&grid_lat[i],4,1,fp);
#ifdef WORDS_BIGENDIAN
#else
                swapFloat(&grid_lat[i]);
#endif

		fseek(fp,offset+4,SEEK_SET);
		fread(&grid_lon[i],4,1,fp);
#ifdef WORDS_BIGENDIAN
#else
                swapFloat(&grid_lon[i]);
#endif

		fseek(fp,offset+8,SEEK_SET);
		fread(&pixel_number_of_centre[i],2,1,fp);
#ifdef WORDS_BIGENDIAN
#else
                swapShortInt(&pixel_number_of_centre[i]);
#endif

		fseek(fp,offset+10,SEEK_SET);
                fread(&line_number_of_centre[i],2,1,fp);
#ifdef WORDS_BIGENDIAN
#else
                swapShortInt(&line_number_of_centre[i]);
#endif

		i++;
	}



	for (j=0;j<DIMEN;j+=Width){
		inverselat=(j-line_number_of_centre[0])*(mercatorProjLat(grid_lat[no_grid_pt-1])- mercatorProjLat(grid_lat[0]))/(line_number_of_centre[no_grid_pt-1]-line_number_of_centre[0]) + mercatorProjLat(grid_lat[0]);
		lat=inverseMercatorProjLat(inverselat);
		for (i=0;i<DIMEN;i+=Width){
			lon=(i-pixel_number_of_centre[0])*(grid_lon[no_grid_pt-1]-grid_lon[0])/
				 (pixel_number_of_centre[no_grid_pt-1]-pixel_number_of_centre[0])+grid_lon[0];  
			fseek(fp,offset_imageData+width_image_pixel*j+i,SEEK_SET);
		    	fread(&pixel_value,1,1,fp);

		  /* NOTE that we used [j][i] instead of [i][j] -  transpose for going back to Fortran */
		  /* when it return, [0][0] should be at sw corner */
			image_data[DIMEN-j-1][i] = pixel_value;
		}
	}

	/* read calibration table */
	for (i=0; i<DIMEN_CAL; i++){
		fseek(fp,offset_cal+4*256+4*i,SEEK_SET);
		fread(&cal_value,4,1,fp);
#ifdef WORDS_BIGENDIAN
#else
		swapFloat(&cal_value);
#endif
		calibration_table[i] = cal_value;
	}

 
	fclose(fp);
/*
	fclose(ft);
*/

	return;
}

float mercatorProjLat(float latitude){
	return log(1/tan((45-latitude/2)*Pi/180));
}

float inverseMercatorProjLat(float pt){
	return 90-(360/Pi)*atan(exp(-pt));
}

#ifdef WORDS_BIGENDIAN
#else
void swapShortInt(short int *n){
   unsigned char *cptr,tmp;
   cptr = (unsigned char *)n;
   tmp = cptr[0];
   cptr[0] = cptr[1];
   cptr[1] =tmp;

}

void swapInt(int *n){
        unsigned char *cprt,tmp;
        cprt=(unsigned char *)n;
        tmp=cprt[0];
        cprt[0]=cprt[3];
        cprt[3]=tmp;
        tmp=cprt[1];
        cprt[1]=cprt[2];
        cprt[2]=tmp;
}

void swapFloat(float *n){
   unsigned char *cptr,tmp;
   cptr = (unsigned char *)n;
   tmp = cptr[0];
   cptr[0] = cptr[3];
   cptr[3] =tmp;
   tmp = cptr[1];
   cptr[1] = cptr[2];
   cptr[2] = tmp;
}
#endif

