#!/bin/sh 

# Annotate a simulated weather image
# This script should run in the same directory as the label files
# Variables $ILOC, $MODE_ALLSKY, $POINT, and $YDISP are inherited from environment

  echo "start annotate_allsky.sh for $ILOC"
  echo "MODE_ALLSKY = $MODE_ALLSKY"
  echo "POINT = $POINT"

  ALT_SCALE=`head -3 label2.$ILOC | tail -1 | cut -c17-23`
  AZI_SCALE=`head -3 label2.$ILOC | tail -1 | cut -c24-30`
  LATLON="`head -1 label2.$ILOC`"
  LATLON=`echo $LATLON | sed 's/^[ \t]*//'` # remove leading spaces

  echo "AZI_SCALE = $AZI_SCALE"

# Annotate Model
  if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
    IMGGEOM=`identify allsky_polar_$ILOC.png | awk '{print tolower($3)}'`
    echo "polar IMGGEOM = $IMGGEOM"
    convert -fill white -annotate +5+20  "Simulated" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
  fi

  if test "$MODE_ALLSKY" = "cyl" || test "$MODE_ALLSKY" = "both"; then
      if test $AZI_SCALE == 0.10; then
          convert -fill yellow -annotate +19+447 "Simulated"  -pointsize 16 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      elif test $AZI_SCALE == 0.20; then
          convert -fill yellow -annotate +19+179 "Simulated"  -pointsize 14 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      elif test $AZI_SCALE == 0.25; then
          convert -fill yellow -annotate +15+$YDISP  "Simulated"  -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      elif test $AZI_SCALE == 0.50; then
#         convert -fill yellow -annotate +15+$YDISP  "Simulated"  -pointsize 12 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
          echo " "
      else
          convert -fill yellow -annotate +15+20  "Simulated"  -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      fi
  fi
  
# Annotate Time and GHI
  ATIME=`head -2 label.$ILOC | tail -1`
  ATIME=`echo $ATIME | sed 's/^[ \t]*//'` # remove leading spaces
  GHIUNITS=W/m^2
  GHI=`head -5 label2.$ILOC | tail -1`$GHIUNITS
  GHI=`echo $GHI | sed 's/^[ \t]*//'`     # remove leading spaces
  echo "Annotate Time $ATIME and GHI $GHI"
  if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
   if test "$IMGGEOM" = "511x511"; then
    POINTP=16
    echo 'convert -fill white -annotate +360+503  "$ATIME" -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png'
          convert -fill white -annotate +360+503  "$ATIME" -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   elif test "$IMGGEOM" = "1023x1023"; then
    POINTP=22
    convert -fill white -annotate +770+1003 "$ATIME" -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   elif test "$IMGGEOM" = "1535x1535"; then
    convert -fill white -annotate +1080+1527 "$ATIME" -pointsize 15 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   elif test "$IMGGEOM" = "2047x2047"; then
    convert -fill white -annotate +1440+2039 "$ATIME" -pointsize 15 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   elif test "$IMGGEOM" = "2559x2559"; then
    convert -fill white -annotate +2300+2531 "$ATIME" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   else
    convert -fill white -annotate +1080+1527 "$ATIME" -pointsize 15 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   fi
  fi

  if test "$MODE_ALLSKY" = "cyl" || test "$MODE_ALLSKY" = "both"; then
    if test $AZI_SCALE == 0.10; then
      convert -fill yellow -annotate +525+447 "$ATIME"  -pointsize 16 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.20; then
      convert -fill yellow -annotate +918+179 "$ATIME" -pointsize 14 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.25; then
      XDISP1=480
      XDISP2=815
      if test "$GHI" != ""; then
        echo  convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI" -annotate +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
              convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI" -annotate +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      else
        echo  convert -fill yellow -annotate  +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
              convert -fill yellow -annotate  +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
        echo  convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI"   -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
              convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI"   -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      fi
    elif test $AZI_SCALE == 0.50; then
#     convert -fill yellow -annotate +1550+20 "$ATIME" -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      if test "$GHI" != ""; then
        XDISP1=900
        XDISP2=1550
        convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI" -annotate +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      else
        echo 'convert -fill yellow -annotate  +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png'
              convert -fill yellow -annotate  +$XDISP2+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
        echo 'convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI"   -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png'
              convert -fill yellow -annotate  +$XDISP1+$YDISP "$GHI"   -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      fi
    elif test "$IMGGEOM" = "5761x2881"; then
      XDISP2=3000
      convert -fill white -annotate +$XDISP2+$YDISP "$ATIME" -pointsize 25 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    else
      convert -fill yellow -annotate +725+20 "$ATIME"  -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    fi
  fi
  
# Annotate Lat/Lon
# Try and strip off leading blanks?
# LATLON="39.99 -105.26"
  echo "Annotate Lat/Lon $LATLON"
  if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
   if test "$IMGGEOM" = "511x511"; then
    echo "convert -fill white -annotate +363+20 "$LATLON" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png"
          convert -fill white -annotate +363+20 "$LATLON" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   elif test "$IMGGEOM" = "1023x1023"; then
    echo "convert -fill white -annotate +875+20 "$LATLON" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png"
          convert -fill white -annotate +875+20 "$LATLON" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   elif test "$IMGGEOM" = "2559x2559"; then
    echo "convert -fill white -annotate +2300+20 "$LATLON" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png"
          convert -fill white -annotate +2300+20 "$LATLON" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   else
    echo "convert -fill white -annotate +1387+20 "$LATLON" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png"
          convert -fill white -annotate +1387+20 "$LATLON" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   fi
  fi

  if test "$MODE_ALLSKY" = "cyl" || test "$MODE_ALLSKY" = "both"; then
    if test $AZI_SCALE == 0.10; then
      convert -fill yellow -annotate +820+447 "$LATLON" -pointsize 16 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.20; then
      convert -fill yellow -annotate +1434+179 "$LATLON" -pointsize 14 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.25; then
      convert -fill yellow -annotate +1227+$YDISP "$LATLON" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.50; then
      convert -fill yellow -annotate +1840+20 "$LATLON" -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    else
      convert -fill yellow -annotate +920+20  "$LATLON" -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    fi
  fi

# Annotate Field
  if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
    if test "$GHI" != ""; then
      LL=$GHI
      LLPOINT=16
    else
      LL=All-sky
      LLPOINT=18
    fi

    if test "$IMGGEOM" = "511x511"; then
      convert -fill white -annotate +20+500  "$LL"   -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    elif test "$IMGGEOM" = "1023x1023"; then
      convert -fill white -annotate +20+1003 "$LL"   -pointsize $POINTP allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    elif test "$IMGGEOM" = "1535x1535"; then
      convert -fill white -annotate +20+1524 "$LL"   -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    elif test "$IMGGEOM" = "2047x2047"; then
      convert -fill white -annotate +20+1524 "$LL"   -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    elif test "$IMGGEOM" = "2559x2559"; then
      convert -fill white -annotate +20+2531 "$LL"   -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    else
      convert -fill white -annotate +20+1524 "$LL"   -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    fi
  fi
