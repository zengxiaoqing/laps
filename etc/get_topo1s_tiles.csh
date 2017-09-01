#!/bin/csh 

# Obtain 1s topo data for 1 degree tiles in domain

# Setup of this script:
# Edit GEOG variable and foreach LAT/LON ranges
# Download GDAL software at http://www.gdal.org/

setenv GEOG /Users/albers/software/geog

setenv DATA1 $GEOG/arcgrid 
setenv DATA2 $GEOG/topo_1s

mkdir -p $DATA1
mkdir -p $DATA2

cd $DATA1

foreach LAT (n40 n41 n42 n43 n44 n45 n46)
#foreach LAT (n42)
foreach LON (w107 w108 w109 w110 w111 w112 w113 w114 w115 w116)
    setenv LATLON $LAT$LON

#   setenv FILE USGS_NED_1_$LATLON\_ArcGrid     
    setenv FILE $LATLON                       

    rm -f $FILE 

    if (! -e $FILE) then
#       wget https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/ArcGrid/$FILE                            
        wget https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/ArcGrid/$FILE.zip
        unzip -o $FILE.zip
    else
        echo $FILE already exists...
    endif

    setenv TOPODIR grd$LATLON\_1

    setenv BASETOPO img$LATLON\_1
    setenv BASETOPO_SCALED $BASETOPO\_SCALED

    if (! -e $DATA1/$TOPODIR/w001001.adf) then
        echo "Error - $DATA1/$TOPODIR/w001001.adf does not exist"
        exit
    else
        echo "Processing $DATA1/$TOPODIR/w001001.adf"
    endif
    
    gdal_translate -of PNM -ot UInt16                   $DATA1/$TOPODIR/w001001.adf   $DATA2/$LATLON.ppm              
    gdal_translate -of PNM -ot Byte   -scale 1400 2000  $DATA1/$TOPODIR/w001001.adf   $DATA2/$LATLON\_SCALED.ppm
    gdal_translate -of AAIGRID                          $DATA1/$TOPODIR/w001001.adf   $DATA2/$LATLON.asc

    convert -compress none $DATA2/$LATLON.ppm         $DATA2/$LATLON.ppm              
    convert -compress none $DATA2/$LATLON\_SCALED.ppm $DATA2/$LATLON\_SCALED.ppm

end
end

echo "end of script - look for results in $DATA2"
