#!/bin/sh 

# This script should run in the same directory as the label files
# Variables $ILOC, $MODE_ALLSKY, $POINT, and $YDISP are inherited from environment

  echo "start annotate_allsky.sh for $ILOC"
  echo "MODE_ALLSKY = $MODE_ALLSKY"
  echo "POINT = $POINT"

  ALT_SCALE=`head -3 label2.$ILOC | tail -1 | cut -c17-23`
  AZI_SCALE=`head -3 label2.$ILOC | tail -1 | cut -c24-30`
  
# Annotate Time and GHI
  ATIME=`head -2 label.$ILOC | tail -1`
  ATIME=`echo $ATIME | sed 's/^[ \t]*//'` # remove leading spaces
  GHIUNITS=W/m^2
  GHI=`head -5 label2.$ILOC | tail -1`$GHIUNITS
  GHI=`echo $GHI | sed 's/^[ \t]*//'`     # remove leading spaces
  echo "Annotate Time $ATIME and GHI $GHI"
  if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
   IMGGEOM=`identify allsky_polar_$ILOC.png | awk '{print tolower($3)}'`
   if test "$IMGGEOM" = "511x511"; then
    echo 'convert -fill white -annotate +360+503  "$ATIME" -pointsize 15 allsky_polar_$ILOC.png allsky_polar_$ILOC.png'
          convert -fill white -annotate +360+503  "$ATIME" -pointsize 15 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   elif test "$IMGGEOM" = "1023x1023"; then
    convert -fill white -annotate +720+1015 "$ATIME" -pointsize 15 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
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
    else
      convert -fill yellow -annotate +725+20 "$ATIME"  -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    fi
  fi
  
