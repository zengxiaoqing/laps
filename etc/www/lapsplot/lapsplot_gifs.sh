#!/bin/sh 

#umask 000

echo "start lapsplot_gifs.sh"

uname -a

proc=$1
WINDOW=$2
LAPS_GIFS=$3
delay=$4
EXE_DIR=$5
export LAPS_DATA_ROOT=$6
export NCARG_ROOT=$7

# Choices are "yes", "no", or the number of images if a montage is desired
animate=$8

RESOLUTION=$9

SCRATCH_DIR=$LAPS_GIFS/scratch
LAPSPLOT_IN=$SCRATCH_DIR/lapsplot.in_$proc

uscore="_"
MACHINE=`uname -s`
NODE=`uname -n`

#export NCARG_ROOT=/usr/local/apps/ncarg-4.2.2-pgi
SUPMAP_DATA_DIR=/home/elvis/mcdonald/data/supmap/
#alias ctrans '/usr/local/apps/ncarg-4.0.1/bin/ctrans  -verbose'

echo "proc ="$proc
echo "WINDOW ="$WINDOW
echo "LAPS_GIFS = "$LAPS_GIFS
echo "SCRATCH_DIR =  $SCRATCH_DIR"
echo "EXE_DIR =  $EXE_DIR"
echo "NCARG_ROOT ="$NCARG_ROOT
echo "LAPS_DATA_ROOT ="$LAPS_DATA_ROOT
echo "latest ="$latest
echo "RESOLUTION ="$RESOLUTION
echo "LAPSPLOT_IN ="$LAPSPLOT_IN
echo "animate ="$animate
echo "delay ="$delay

echo " "
echo "setting ulimit"
ulimit -t 1000
ulimit -t

#mkdir -p /scratch/lapb/www
mkdir -p $SCRATCH_DIR/$proc
cd $SCRATCH_DIR/$proc

#EXE_DIR=/usr/nfs/lapb/parallel/laps/bin

date -u
echo "Running $EXE_DIR/lapsplot.exe < $LAPSPLOT_IN on $MACHINE $NODE"
rm -f gmeta
#$EXE_DIR/lapsplot.exe                                          < $LAPS_GIFS/lapsplot.in
$EXE_DIR/lapsplot.exe                                           < $LAPSPLOT_IN

if test "$NCARG_ROOT" = "allsky"; then 
 echo "allsky option"
 echo "$SCRATCH_DIR/$proc"
 ls -l $SCRATCH_DIR/$proc

 RESOLUTION_POLAR=""
 RESOLUTION_CYL=""

 if test "$RESOLUTION" = "180p" || test "$RESOLUTION" = "180pr"; then
      MODE_ALLSKY=polar
      RESOLUTION_POLAR=$RESOLUTION
 elif test "$RESOLUTION" = "180c"; then
      MODE_ALLSKY=cyl   
      RESOLUTION_CYL=$RESOLUTION
 elif test "$RESOLUTION" = "360c"; then
      MODE_ALLSKY=cyl   
      RESOLUTION_CYL=$RESOLUTION
 elif test "$RESOLUTION" = "360p" || test "$RESOLUTION" = "360pr"; then
      MODE_ALLSKY=polar 
      RESOLUTION_POLAR=$RESOLUTION
 elif test "$RESOLUTION" = "360b"; then
      MODE_ALLSKY=both   
      RESOLUTION_POLAR=360p
      RESOLUTION_CYL=360c
 elif test "$RESOLUTION" = "180b"; then
      MODE_ALLSKY=both   
      RESOLUTION_POLAR=180p
      RESOLUTION_CYL=180c
 elif test "$RESOLUTION" = "360p180c"; then
      MODE_ALLSKY=both   
      RESOLUTION_POLAR=360p
      RESOLUTION_CYL=180c
 else
      RESOLUTION=360p
      MODE_ALLSKY=polar
 fi

 NLOC=`head -5      $LAPSPLOT_IN | tail -1`
 echo "NLOC from input file $LAPSPLOT_IN is $NLOC"

#NLOC=1

 for (( i=1; i<=$NLOC; i++ )); do
  export ILOC=00$i
  echo "Welcome $ILOC times"

  echo "MODE_ALLSKY = $MODE_ALLSKY    RESOLUTION = $RESOLUTION"
  echo "RESOLUTION_CYL = $RESOLUTION_CYL    RESOLUTION_POLAR = $RESOLUTION_POLAR"

  MINALT=`head -2 label2.$ILOC | tail -1 | cut -c1-8`
  MAXALT=`head -2 label2.$ILOC | tail -1 | cut -c9-16`
  MINAZI=`head -2 label2.$ILOC | tail -1 | cut -c17-24`
  MAXAZI=`head -2 label2.$ILOC | tail -1 | cut -c25-32`
  NI_CYL=`head -2 label2.$ILOC | tail -1 | cut -c33-40`
  NJ_CYL=`head -2 label2.$ILOC | tail -1 | cut -c41-48`
  ALT_SCALE=`head -3 label2.$ILOC | tail -1 | cut -c17-23`
  AZI_SCALE=`head -3 label2.$ILOC | tail -1 | cut -c24-30`
  LATLON="`head -1 label2.$ILOC`"
  LATLON=`echo $LATLON | sed 's/^[ \t]*//'` # remove leading spaces

  echo "MINALT/MAXALT/MINAZI/MAXAZI/NI_CYL/NJ_CYL = $MINALT $MAXALT $MINAZI $MAXAZI $NI_CYL $NJ_CYL"
  echo "ALT_SCALE/AZI_SCALE = $ALT_SCALE,$AZI_SCALE"
  echo "LATLON=$LATLON"

  if test "$LATLON" = "38.84 -105.04"; then
    echo "Site is Pikes"
    YDISP=135
    POINT=12
  else
    echo "Site is not Pikes"
    YDISP=20
    POINT=16
  fi

  echo "YDISP=$YDISP"

  if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
    echo " "
#   echo "Will run IDL polar conversion to PNG"
#   rm -f allsky*.pro; ln -s /home/fab/albers/ast/skyglow/allsky.pro allsky.pro
#   echo allsky | /usr/local/share/rsi/idl/bin/idl
    echo "Will run convert -flip allsky_rgb_polar_$ILOC.ppm allsky_polar_$ILOC.png"
    convert -flip allsky_rgb_polar_$ILOC.ppm allsky_polar_$ILOC.png
  fi

  if test "$MODE_ALLSKY" = "cyl" || test "$MODE_ALLSKY" = "both"; then
    WINDOW=$NI_CYL # Resetting WINDOW
#   export ALLSKY_JDIM=$WINDOW
    export ALLSKY_IDIM=$NJ_CYL
    export ALLSKY_JDIM=$NI_CYL
    echo " "
#   echo "Will run IDL cyl conversion to PNG: ALLSKY_IDIM/ALLSKY_JDIM is $ALLSKY_IDIM $ALLSKY_JDIM"
#   rm -f allsky*.pro; ln -s /home/fab/albers/ast/skyglow/allsky_cyl.pro allsky_cyl.pro
#   echo allsky_cyl | /usr/local/share/rsi/idl/bin/idl
    echo "Will run convert -flip allsky_rgb_cyl_$ILOC.ppm allsky_cyl_$ILOC.png"
    convert -flip allsky_rgb_cyl_$ILOC.ppm allsky_cyl_$ILOC.png
    if test $NJ_CYL -ge 1441; then
#   if test $NJ_CYL -ge  750; then
      echo "cylindrical image is of sufficient size - no enlargement"
    elif test $NJ_CYL -ge 750; then  
      echo "enlarge cylindrical image by 150%"
      convert -resize 150% allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    else
      echo "enlarge cylindrical image by 300%"
      convert -resize 300% allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    fi
  fi

  echo " "

  DIRCYL=South
  DIRCYLL=East
  DIRCYLR=West

# Other orientations
  if test "$RESOLUTION_POLAR" = "180pr"; then
      convert allsky_polar_$ILOC.png -rotate 180 -flop allsky_polar_$ILOC.png
      DIR1=SE
      DIR2=NE
      DIR3=NW
      DIR4=SW
  fi
  if test "$RESOLUTION_POLAR" = "360pr"; then # default           
      convert allsky_polar_$ILOC.png -flop allsky_polar_$ILOC.png
      DIR1=NW
      DIR2=SW
      DIR3=SE
      DIR4=NE
  fi
  if test "$RESOLUTION_POLAR" = "360p"; then # flip left/right
      DIR1=NE
      DIR2=SE
      DIR3=SW
      DIR4=NW
  fi
  if test "$RESOLUTION_POLAR" = "180p"; then # rotate and flip left/right
      convert allsky_polar_$ILOC.png -rotate 180 allsky_polar_$ILOC.png
#     convert allsky_polar_$ILOC.png -flop allsky_polar_$ILOC.png
      DIR1=SW
      DIR2=NW
      DIR3=NE
      DIR4=SE
  fi
  if test "$RESOLUTION_CYL" = "360c"; then # roll horizontally by half the image
      if test $AZI_SCALE == 0.25; then
          convert allsky_cyl_$ILOC.png -roll  +720+0 allsky_cyl_$ILOC.png
      elif test $AZI_SCALE == 0.50; then
          convert allsky_cyl_$ILOC.png -roll +1080+0 allsky_cyl_$ILOC.png
      else
          convert allsky_cyl_$ILOC.png -roll  +540+0 allsky_cyl_$ILOC.png
      fi
      DIRCYL=North
      DIRCYLL=West
      DIRCYLR=East
  fi

# Annotate Model
  if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
    IMGGEOM=`identify allsky_polar_$ILOC.png | awk '{print tolower($3)}'`
    echo "polar IMGGEOM = $IMGGEOM"
    convert -fill white -annotate +5+20  "LAPS Simulated" -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
  fi

  if test "$MODE_ALLSKY" = "cyl" || test "$MODE_ALLSKY" = "both"; then
      if test $AZI_SCALE == 0.10; then
          convert -fill yellow -annotate +19+447 "LAPS Simulated"  -pointsize 16 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      elif test $AZI_SCALE == 0.20; then
          convert -fill yellow -annotate +19+179 "LAPS Simulated"  -pointsize 14 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      elif test $AZI_SCALE == 0.25; then
          convert -fill yellow -annotate +15+$YDISP  "LAPS Simulated"  -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      else
          convert -fill yellow -annotate +15+20  "LAPS Simulated"  -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      fi
  fi

# Annotate Time
  ATIME=`head -2 label.$ILOC | tail -1`
  ATIME=`echo $ATIME | sed 's/^[ \t]*//'` # remove leading spaces
  echo "Annotate Time $ATIME"
  if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
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
      convert -fill yellow -annotate  +815+$YDISP "$ATIME" -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.50; then
      convert -fill yellow -annotate +1550+20 "$ATIME" -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
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
    if test "$IMGGEOM" = "511x511"; then
      convert -fill white -annotate +20+500  "All-Sky"   -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    elif test "$IMGGEOM" = "1023x1023"; then
      convert -fill white -annotate +20+1012 "All-Sky"   -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    else
      convert -fill white -annotate +20+1524 "All-Sky"   -pointsize 18 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    fi
  fi

# if test "$MODE_ALLSKY" = "cyl" || test "$MODE_ALLSKY" = "both"; then
#   convert -fill white -annotate +20+500 "All-Sky"   -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
# fi

# Annotate Directions
  echo "Annotate Directions"
  if test "$MODE_ALLSKY" = "polar" || test "$MODE_ALLSKY" = "both"; then
   if test "$IMGGEOM" = "511x511"; then
    convert -fill white -annotate +55+60     "$DIR1"            -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    convert -fill white -annotate +40+450    "$DIR2"            -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    convert -fill white -annotate +440+450   "$DIR3"            -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
    convert -fill white -annotate +435+60    "$DIR4"            -pointsize 20 allsky_polar_$ILOC.png allsky_polar_$ILOC.png
   fi
  fi

  if test "$MODE_ALLSKY" = "cyl" || test "$MODE_ALLSKY" = "both"; then
    if test $AZI_SCALE == 0.10; then
      convert -fill orange -annotate  +1050+447   "East"          -pointsize 16 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
    elif test $AZI_SCALE == 0.20; then
      convert -fill orange -annotate  +776+179  "$DIRCYL"          -pointsize 14 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      echo "WINDOW is larger "$WINDOW
    elif test $AZI_SCALE == 0.25; then
      convert -fill orange -annotate  +693+$YDISP   "$DIRCYL"          -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      convert -fill orange -annotate  +339+$YDISP   "$DIRCYLL"         -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      convert -fill orange -annotate  +1055+$YDISP  "$DIRCYLR"         -pointsize $POINT allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      echo "WINDOW is large "$WINDOW
    elif test $AZI_SCALE == 0.50; then
      convert -fill orange -annotate +1040+20   "$DIRCYL"          -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      echo "WINDOW is large "$WINDOW
    else
      convert -fill orange -annotate +520+20    "$DIRCYL"          -pointsize 20 allsky_cyl_$ILOC.png allsky_cyl_$ILOC.png
      echo "WINDOW is small "$WINDOW
    fi
  fi

# convert allsky_polar_ILOC.png $SCRATCH_DIR/gmeta_$proc.$ILOC.gif
  if test "$MODE_ALLSKY" = "cyl"; then
      cp allsky_cyl_$ILOC.png $SCRATCH_DIR/gmeta_$proc.$ILOC.png
      cp allsky_cyl_$ILOC.png $SCRATCH_DIR/gmeta_$proc.png # for on-the-fly page
      ls -l $SCRATCH_DIR/gmeta_$proc.$ILOC.png
  fi
  if test "$MODE_ALLSKY" = "polar"; then
      cp allsky_polar_$ILOC.png $SCRATCH_DIR/gmeta_$proc.$ILOC.png
      cp allsky_polar_$ILOC.png $SCRATCH_DIR/gmeta_$proc.png # for on-the-fly page
      ls -l $SCRATCH_DIR/gmeta_$proc.$ILOC.png
  fi

# handle "both" option?

  if test "$MODE_ALLSKY" = "both"; then
      cp allsky_cyl_$ILOC.png $SCRATCH_DIR/gmeta_$proc.$ILOC.cyl.png
      cp allsky_polar_$ILOC.png $SCRATCH_DIR/gmeta_$proc.$ILOC.polar.png
      ls -l $SCRATCH_DIR/gmeta_$proc.$ILOC.cyl.png
      ls -l $SCRATCH_DIR/gmeta_$proc.$ILOC.polar.png
  fi

  ls -l $SCRATCH_DIR/$proc

  ext3=png

 done # ILOC loop

else
  echo "regular ncar graphics option"
  pwd
  ls -l $SCRATCH_DIR/$proc/gmeta
  if test ! -e gmeta; then
    echo "ERROR: gmeta does not exist in lapsplot_gifs.sh - exiting"
    exit
  fi

  if test "$MACHINE" = "AIX"; then

#   Combination for IBM
    ext1=avs
    ext2=x
    ext3=gif
    netpbm=no

#   CTRANS=/usr/local/apps/ncarg-4.0.1/bin/ctrans

  else

    netpbm=yes

#   Best combination for LINUX
    ext1=avs
    ext2=avs
    ext3=gif

#   The ones below will run but produce fewer colors in color images for LINUX

#   ext1=sun
#   ext2=sun

#   ext1=xwd
#   ext2=xwd

#   Note that sgi will not work in LINUX since we are using gmeta files with WINDOW/RESOLUTION set
#   ext1=sgi
#   ext2=sgi

#   ext1=sun
#   ext2=gif

#   CTRANS=/usr/local/apps/ncarg-4.2.2-pgi/bin/ctrans

  fi

  CTRANS=$NCARG_ROOT/bin/ctrans


# /usr/local/apps/ncarg-4.0.1/bin/ctrans -verbose -d avs -window $WINDOW -resolution $RESOLUTION gmeta > $SCRATCH_DIR/gmeta_$proc.x
# ctrans -d avs -window 0.0:0.08:1.0:0.92 -resolution 610x512 gmeta > $SCRATCH_DIR/gmeta_$proc.x
# /usr/local/apps/ncarg-4.0.1/bin/ctrans -verbose -d $ext1 -window $WINDOW -resolution $RESOLUTION gmeta > $SCRATCH_DIR/gmeta_temp_$proc.$ext2
# $NCARG_ROOT/bin/ctrans -verbose -d $ext1 -window $WINDOW -resolution $RESOLUTION $SCRATCH_DIR/gmeta > $SCRATCH_DIR/gmeta_temp_$proc.$ext2

# $CTRANS -verbose -d $ext1 -window $WINDOW -resolution $RESOLUTION $SCRATCH_DIR/$proc/gmeta > $SCRATCH_DIR/gmeta_temp_$proc.$ext2

# ls -l $SCRATCH_DIR/gmeta_temp_$proc.$ext2

  date -u

  echo "lapsplot_gifs.sh: netpbm = $netpbm"

# numimages=`ls -1 *.gif | wc -l`
# echo "numimages = $numimages"

# We assume we are running this script in LINUX and convert will not properly do AVS X on LINUX
  if test "$netpbm" = "yes" && test "$animate" = "no"; then 
    date
#   echo "Running $NCARG_ROOT/bin/ctrans | netpbm to make gmeta_$proc.gif file"
    echo "Running $NCARG_ROOT/bin/ctrans -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta | rasttopnm | ppmtogif > $SCRATCH_DIR/gmeta_$proc.gif"
    $NCARG_ROOT/bin/ctrans -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta | rasttopnm | ppmtogif > $SCRATCH_DIR/gmeta_$proc.gif

    date -u

#   Cleanup
    echo "Cleanup"
    mv gmeta $SCRATCH_DIR/gmeta_$proc.gm;  cd ..; rmdir $SCRATCH_DIR/$proc &

  elif test "$netpbm" = "yes" && test "$animate" != "no"; then 
    date
    echo "Running $NCARG_ROOT/bin/ctrans -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta > $SCRATCH_DIR/$proc/gmeta_$proc.sun"
    $NCARG_ROOT/bin/ctrans -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta > $SCRATCH_DIR/$proc/gmeta_$proc.sun

#   Convert multiframe raster image to animated gif
    $NCARG_ROOT/bin/rassplit gmeta_$proc.sun

#   Convert sun to gif images so convert works better on new server
    for file in `ls gmeta_$proc.*.sun`; do
        ls $file
        rasttopnm $file | ppmtogif > $file.gif
    done

    ls -l gmeta_$proc.*.sun.gif

#   Make animation or montage
    if test "$animate" != "no"; then
      numimages=`ls -1 *.$ext3 | wc -l`
      echo "numimages = $numimages"

      echo " "
      echo "Listing of $SCRATCH_DIR/$proc animation images"
      ls -1r gmeta*.$ext3 | tee files.txt

      if test "$animate" = "yes"; then
        echo "convert -delay $delay -loop 0 *.$ext3                 $SCRATCH_DIR/gmeta_$proc.gif"
              convert -delay $delay -loop 0 *.$ext3 $file.$ext3     $SCRATCH_DIR/gmeta_$proc.gif

      else # make montage instead of animation, $animate is the number of images
#       numimages=`ls -1 *.gif | wc -l`
#       echo "numimages = $numimages"

        nmontage=$animate

        echo "nmontage = $nmontage"

        montage_file=$SCRATCH_DIR/montage_$proc.sh

        if test -r "$montage_file"; then

          echo "running montage file: $montage_file"
          cat $montage_file
          /bin/sh $montage_file
          rm -f $montage_file

          echo " "
          echo "Listing of $SCRATCH_DIR/$proc animation images"
          ls -1r gmeta_*_*.$ext3 | tee files.txt
          rm -f *sun*.$ext3

        else

          x20=x20
          x=x

          if test "$numimages" == "3"; then # single row
            echo "making single row"
            echo "montage *.$ext3 -mode Concatenate -tile $nmontage$x20 $SCRATCH_DIR/gmeta_$proc.gif"
                  montage *.$ext3 -mode Concatenate -tile $nmontage$x20 $SCRATCH_DIR/gmeta_$proc.gif
          elif test "$numimages" == "4"; then # double row
            echo "making double row"
            echo "montage *.$ext3 -mode Concatenate -tile 2x2           $SCRATCH_DIR/gmeta_$proc.gif"
                  montage *.$ext3 -mode Concatenate -tile 2x2           $SCRATCH_DIR/gmeta_$proc.gif
          else                              # automatic settings
            echo "making $nmontage (nmontage) columns"
            echo "montage *.$ext3 -mode Concatenate -tile $nmontage$x     $SCRATCH_DIR/gmeta_$proc.gif"
                  montage *.$ext3 -mode Concatenate -tile $nmontage$x     $SCRATCH_DIR/gmeta_$proc.gif
          fi    

        fi

      fi # decide between animation and montage

    fi # make animation or montage (animation != "no")

#   This option may be more direct though it isn't working on the new server
#   echo "convert -delay $delay -loop 0 gmeta_$proc.*.sun $SCRATCH_DIR/gmeta_$proc.gif"
#   convert -delay $delay -loop 0 gmeta_$proc.*.sun $SCRATCH_DIR/gmeta_$proc.gif

    ln -s -f /w3/lapb/looper/files.cgi files.cgi

    date -u

#   Cleanup
    echo "Cleanup"
#   mv gmeta $SCRATCH_DIR/gmeta_$proc.gm;  cd ..; rm -f $SCRATCH_DIR/$proc/gmeta*; rmdir $SCRATCH_DIR/$proc &
    mv gmeta $SCRATCH_DIR/gmeta_$proc.gm;  cd ..; rm -f gmeta; rm -f $SCRATCH_DIR/$proc/*.sun & 

#   echo " "
#   echo "Cleaned up listing of $SCRATCH_DIR/$proc"
#   ls -l $SCRATCH_DIR/$proc

  fi # netpbm / animate test                        

fi # allsky option

chmod 666 $SCRATCH_DIR/gmeta_$proc.$ext3
echo " "
date -u
echo " "
echo "finished lapsplot_gifs.sh"
echo " "

