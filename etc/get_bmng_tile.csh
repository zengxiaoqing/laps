#!/bin/csh 

# Run this after the localization "window" script

#May

#setenv LAPS_DATA_ROOT /data/lapb/operational/laps/data
#setenv LAPS_DATA_ROOT /data/fab/dlaps/projects/roc/hires2
setenv LAPS_DATA_ROOT /Users/albers/data/domains/wrnus_hi
#setenv LAPS_DATA_ROOT /w3/jet/fab/laps_conus
#setenv LAPS_DATA_ROOT /data/fab/dlaps/wrnus
#setenv LAPS_DATA_ROOT /data/fab/projects/muri
#setenv LAPS_DATA_ROOT $1                       
setenv STATIC $LAPS_DATA_ROOT/static
setenv GEOG /Users/albers/data

echo " "
echo "Running command_tile for topo_1s files..."
echo "LAPS_DATA_ROOT is $LAPS_DATA_ROOT"
echo "LATLON domain bounds from llbounds.dat: `cat $STATIC/llbounds.dat`"

cd $GEOG/nasa

# Obtain and process a single Blue Marble tile and make a crop for the domain
# Image section is 90x90 degrees at 240 pixels per degree
# Overall browse page: http://visibleearth.nasa.gov/view_cat.php?categoryID=1484
# Blue Marble color: http://earthobservatory.nasa.gov/Features/BlueMarble/bmng.pdf
#
# Original values for 500m domain
# shift of 17000 = 70.8333 deg
# shift of 11000 = 45.8333 deg
# width of 2000  =  8.3333 deg
# lower lat      = 35.8333 
# upper lat      = 44.1666
# lower lon      = -109.166667
# upper lon      = -100.833333

# Set new values for domain of interest
# Dynamic values have a perimeter of 0.2 around 'llbounds.dat' values
# Longitude increases to the east
setenv PERIMETER 0.2
setenv ULAT `head -1 $STATIC/llbounds.dat | tail -1`
setenv ULAT `echo $ULAT $PERIMETER | awk '{print ($1+$2)}'`
setenv LLAT `head -2 $STATIC/llbounds.dat | tail -1`
setenv LLAT `echo $LLAT $PERIMETER | awk '{print ($1-$2)}'`
setenv ULON `head -3 $STATIC/llbounds.dat | tail -1`
setenv ULON `echo $ULON $PERIMETER | awk '{print ($1+$2)}'`
setenv LLON `head -4 $STATIC/llbounds.dat | tail -1`
setenv LLON `echo $LLON $PERIMETER | awk '{print ($1-$2)}'`

echo "PERIMETER is $PERIMETER"
echo "LAT crop bounds (north to south) are $ULAT $LLAT"
echo "LON crop bounds (west to east)   are $LLON $ULON"

setenv TILE A1  # A1 North America, B1 North Atlantic, C1 Europe, D1 Asia 
foreach MN (01 02 03 04 05 06 07 08 09 10 11 12)

  if ($MN == 01) then
    setenv ID 73938
    setenv IDT 73000
  else if ($MN == 02) then
    setenv ID 73967
    setenv IDT 73000
  else if ($MN == 03) then
    setenv ID 73992
    setenv IDT 73000
  else if ($MN == 04) then
    setenv ID 74017
    setenv IDT 74000
  else if ($MN == 05) then
    setenv ID 74042
    setenv IDT 74000
  else if ($MN == 06) then
    setenv ID 76487
    setenv IDT 76000
  else if ($MN == 07) then
    setenv ID 74092
    setenv IDT 74000
  else if ($MN == 08) then
    setenv ID 74117
    setenv IDT 74000
  else if ($MN == 09) then
    setenv ID 74142
    setenv IDT 74000
  else if ($MN == 10) then
    setenv ID 74167
    setenv IDT 74000
  else if ($MN == 11) then
    setenv ID 74192
    setenv IDT 74000
  else if ($MN == 12) then
    setenv ID 74218
    setenv IDT 74000
  endif

  setenv FILENAME world.2004$MN.3x21600x21600.$TILE.png

  if ($TILE == A1) then   # North America
    setenv WLONTILE -180.
  else                    # Asia
    setenv WLONTILE   90.
  endif

  echo "West Longitude of Tile is $WLONTILE"

# echo '"$LLON < $WLONTILE"' > /tmp/compare
  setenv COMPARE_RESULT `echo "$LLON < $WLONTILE" | bc`
  if ($COMPARE_RESULT == 1) then
    echo "west edge is left of the tile - reset LLON"
    setenv LLON $WLONTILE
  endif
  echo "COMPARE_RESULT / LLON is $COMPARE_RESULT / $LLON"

  setenv ICROP `echo $ULON $LLON     | awk '{print     int((($1-$2)*240.)+0.5)}'`
  setenv JCROP `echo $ULAT $LLAT     | awk '{print     int((($1-$2)*240.)+0.5)}'`
  setenv IOFF  `echo $LLON $WLONTILE | awk '{print     int((($1-$2)*240.)+0.5)}'`
  setenv JOFF  `echo $ULAT           | awk '{print     int(((90.-$1 )*240.)+0.5)}'`

  echo "ICROP = $ICROP"
  echo "JCROP = $JCROP"
  echo "IOFF  = $IOFF"
  echo "JOFF  = $JOFF"

  setenv STRING $ICROP\x$JCROP+$IOFF+$JOFF
  echo "STRING = $STRING"

  if (! -e $FILENAME) then
    echo "$FILENAME will be downloaded"
    echo "wget http://eoimages.gsfc.nasa.gov/images/imagerecords/$IDT/$ID/$FILENAME"
          wget http://eoimages.gsfc.nasa.gov/images/imagerecords/$IDT/$ID/$FILENAME          
  else
    echo "$FILENAME already exists"
  endif
  ls -l $FILENAME                               

  setenv CROPDIR $STATIC
  setenv CROPFILE $CROPDIR/world.2004$MN.3x21600x21600.crop.ppm

  echo "convert -crop $STRING -compress none world.2004$MN.3x21600x21600.$TILE.png $CROPFILE"
        convert -crop $STRING -compress none world.2004$MN.3x21600x21600.$TILE.png $CROPFILE

  ls -l $CROPFILE

  convert -rotate 0 $CROPFILE $CROPDIR/bm.crop.png
  ls -l $CROPDIR/bm.crop.png

end

echo "command_tile completion..."
