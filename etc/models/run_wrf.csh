#!/bin/csh

# Conditions:
# METGRID.TBL namelist.wps.bdy namelist.wps.init namelist.input.bdy namelist.input.init, parame.in
# are available at the current directory
# NOTE: METGRID.TBL either directly from WRFSRC or modified by yourself for VVEL
#       namelist.wps.bdy:  modified to fit your domain. 
#       namelist.wps.init: modified to fit your domain.
#                    Important:  Have to specify the TBL paths;
#       namelist.input.bdy: modified for your boundary files;
#       namelist.input.init: modified for using LAPS as initial condition
#
# WRF: WPS, WRFV3 and WRFDA are sucessfully compiled.

# Usage:
# run_wrf.csh lapsdir wrfdir wrfsrc fcsttime
if ($1 == '' || $2 == '' || $3 == '' || $4 == '') then
  echo "Usage: run_wrf.csh lapsdir wrfdir wrfsrc fcsttime"
  exit
endif

# User needs modification of these paths:
setenv BDYFILES "/data/fab/cwb/windsor/archive/data/grids/ruc/iso/grib1/081431"

setenv LAPS_DATA_ROOT $1
setenv WRFDIR $2
setenv WRFSRC $3
setenv FCSTIME $4
setenv YEAR `echo $4 | cut -c1-4`
setenv MONTH `echo $4 | cut -c5-6`
setenv DAY `echo $4 | cut -c7-8`
setenv HOUR `echo $4 | cut -c9-10`

setenv CDIR `pwd`

echo $LAPS_DATA_ROOT $WRFDIR $CDIR

# working directory:
mkdir -p $WRFDIR
cd $WRFDIR
mkdir $FCSTIME
mkdir $FCSTIME/static
mkdir $FCSTIME/wpsprd
mkdir $FCSTIME/wrfprd

# setup static:
/bin/rm -f $FCSTIME/static/GEOGRID.TBL
ln -s $WRFSRC/WPS/geogrid/GEOGRID.TBL.ARW $FCSTIME/static/GEOGRID.TBL
# Copy METGRID.TBL for modification:
echo $CDIR/METGRID.TBL.bdy
/bin/rm -f $FCSTIME/static/METGRID.TBL
cp $CDIR/METGRID.TBL.bdy $FCSTIME/static/METGRID.TBL
/bin/rm -f $FCSTIME/static/lfmpost.nl
ln -s $LAPS_DATA_ROOT/static/lfmpost.nl $FCSTIME/static/lfmpost.nl

# wps:
# run geogrid:
/bin/rm -f $FCSTIME/wpsprd/namelist.wps
cp $CDIR/namelist.wps.bdy $FCSTIME/wpsprd/namelist.wps
cd $FCSTIME/wpsprd
echo "+---------------------+"
echo "|  Start geogrid.exe  |"
echo "+---------------------+"
$WRFSRC/WPS/geogrid.exe
echo "+---------------------+"
echo "| End of geogrid.exe  |"
echo "+---------------------+"
echo ""

# run ungrid:
/bin/rm -f Vtable
ln -s $WRFSRC/WPS/ungrib/Variable_Tables/Vtable.RUCp Vtable
echo $BDYFILES
/bin/rm -f GRIBFILE*
$WRFSRC/WPS/link_grib.csh "$BDYFILES*"
echo "+---------------------+"
echo "|  Start ungrib.exe   |"
echo "+---------------------+"
$WRFSRC/WPS/ungrib.exe >& ungrib.output
echo "+---------------------+"
echo "|  End of ungrib.exe  |"
echo "+---------------------+"
echo ""

# run metgrid:
echo "+----------------------------+"
echo "| Start metgrid.exe for bdy  |"
echo "+----------------------------+"
$WRFSRC/WPS/metgrid.exe 
echo "+----------------------------+"
echo "| End of metgrid.exe for bdy |"
echo "+----------------------------+"
echo ""

# run real.exe for boundary files:
cd ../wrfprd
mv ../wpsprd/met_em* .
/bin/rm -f namelist.input
cp $CDIR/namelist.input.bdy namelist.input
echo "+--------------------------+"
echo "|  Start real.exe for bdy  |"
echo "+--------------------------+"
$WRFSRC/WRFV3/run/real.exe
echo "+--------------------------+"
echo "|   End real.exe for bdy   |"
echo "+--------------------------+"
echo ""
mv wrfbdy_d01 wrfbdy_d01.save
mv wrfinput_d01 wrfinput_d01.save

# LAPS/STMAS initial condition:
cd ../wpsprd
/bin/rm -f "LAPS:$YEAR""-$MONTH""-$DAY""_$HOUR"
ln -s "$LAPS_DATA_ROOT/lapsprd/lapsprep/wps/LAPS:$YEAR""-$MONTH""-$DAY""_$HOUR" .
/bin/rm -f namelist.wps
cp $CDIR/namelist.wps.init namelist.wps
ls $CDIR/METGRID.TBL.vvel
/bin/rm -f ../static/METGRID.TBL
cp $CDIR/METGRID.TBL.init ../static/METGRID.TBL
echo "+------------------------------+"
echo "|  Start metgrid.exe for init  |"
echo "+------------------------------+"
$WRFSRC/WPS/metgrid.exe
echo "+------------------------------+"
echo "| End of metgrid.exe for in-it |"
echo "+------------------------------+"
echo ""
setenv NEWFILE `ls met_em*`

# run real.exe for initial file:
cd ../wrfprd
echo $NEWFILE
/bin/rm -f $NEWFILE
mv ../wpsprd/met_em* .
/bin/rm -f namelist.input
cp $CDIR/namelist.input.init namelist.input
echo "+---------------------------+"
echo "|  Start real.exe for init  |"
echo "+---------------------------+"
$WRFSRC/WRFV3/run/real.exe
echo "+---------------------------+"
echo "|   End real.exe for init   |"
echo "+---------------------------+"
echo ""

# Update boundary conditions: !!! require WRFDA installed
mv wrfinput_d01 wrfvar_output
mv wrfbdy_d01.save wrfbdy_d01
cp $CDIR/parame.in .
echo "+-------------------------+"
echo "|  Start boundary update  |"
echo "+-------------------------+"
$WRFSRC/WRFDA/var/build/da_update_bc.exe
echo "+-------------------------+"
echo "|   End boundary update   |"
echo "+-------------------------+"
echo ""
mv wrfvar_output wrfinput_d01

# run forecasts:
/bin/rm -f LANDUSE.TBL RRTA_DATA
ln -s $WRFSRC/WRFV3/run/LANDUSE.TBL .
ln -s $WRFSRC/WRFV3/run/RRTM_DATA .
echo "+-------------------------+"
echo "|  Start WRF forecasts... |"
echo "+-------------------------+"
$WRFSRC/WRFV3/run/wrf.exe | tee wrfcst.output
cp $CDIR/parame.in .
echo "+-------------------------+"
echo "| End of WRF forecasts... |"
echo "+-------------------------+"
echo ""

echo "*******************************"
echo "*  Forecast script completes  *"
echo "*******************************"


