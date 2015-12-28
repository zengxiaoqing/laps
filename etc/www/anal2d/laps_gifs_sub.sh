#!/bin/sh

#laps_gifs_sub.sh

umask 000

prod=$1
WINDOW=$2
LAPS_ETC=$3
WWW_DIR=$4
LAPS_A9TIME=$5
LAPS_DATA_ROOT=$6
latest=$7
datetime=$8
RESOLUTION=$9

utc_hhmm=`echo $LAPS_A9TIME | cut -c6-9`

SCRATCH_DIR=$WWW_DIR/anal2d

uscore="_"
MACHINE=`uname -n`
PLATFORM=`uname -s`

#export NCARG_ROOT=/usr/local/apps/ncarg-4.0.1
#setenv NCARG_ROOT `cat /usr/nfs/lapb/bin/ncarg_root`

echo "Start laps_gifs_sub.sh..."
echo "PATH=$PATH"

cd
pwd
ls -l etc/ncarg_root
if test -r etc/ncarg_root; then
    echo "Getting NCARG_ROOT from ~/etc/ncarg_root"
    ls -l etc/ncarg_root
    export NCARG_ROOT=`cat etc/ncarg_root`
else
    echo "unable to find ~/etc/ncarg_root file"
fi

#alias ctrans '/usr/local/apps/ncarg-4.0.1/bin/ctrans  -verbose'

echo "prod ="$prod
echo "WINDOW ="$WINDOW
echo "NCARG_ROOT ="$NCARG_ROOT
echo "LAPS_A9TIME ="$LAPS_A9TIME
echo "utc_hhmm ="$utc_hhmm
echo "LAPS_DATA_ROOT ="$LAPS_DATA_ROOT
echo "latest ="$latest
echo "RESOLUTION ="$RESOLUTION
echo "LAPS_ETC = "$LAPS_ETC
echo "WWW_DIR ="$WWW_DIR
echo "SCRATCH_DIR ="$SCRATCH_DIR

cd $WWW_DIR/anal2d

ulimit -t 300

if test -r $EXE_DIR/lapsplot.exe; then
    echo "Running $EXE_DIR/lapsplot.exe"; date -u
else
    echo "ERROR: $EXE_DIR/lapsplot.exe not found"
    exit
fi

#head -2 $LAPS_DATA_ROOT/time/systime.dat | tail -1 | cut -c2-10     > $SCRATCH_DIR/lapsplot.$prod.tmp
echo $LAPS_A9TIME                                                    > $SCRATCH_DIR/lapsplot.$prod.tmp

if test -r $LAPS_DATA_ROOT/static/www/lapsplot.$prod; then
    echo "Input to lapsplot.exe (dataroot) = $LAPS_DATA_ROOT/static/www/lapsplot.$prod"
    cat     $LAPS_DATA_ROOT/static/www/lapsplot.$prod              >> $SCRATCH_DIR/lapsplot.$prod.tmp
elif test -r $SCRATCH_DIR/lapsplot.$prod; then
    echo "Input to lapsplot.exe (local) = $SCRATCH_DIR/lapsplot.$prod"
    cat     $SCRATCH_DIR/lapsplot.$prod                            >> $SCRATCH_DIR/lapsplot.$prod.tmp
else
    echo "ERROR: lapsplot.$prod not found under $LAPS_DATA_ROOT or $SCRATCH_DIR"
    exit
fi

$EXE_DIR/lapsplot.exe                                           < $SCRATCH_DIR/lapsplot.$prod.tmp

pwd; ls -l gmeta

# Test whether we are on NOAA/ESRL/GSD/FAB machines
if test -r /data/fab; then
    netpbm=no
else
    netpbm=yes
fi

# This if block is used to further assess whether we have access to netpbm
if test "$MACHINE" = "headnode.fsl.noaa.gov"; then
    ctransarg=avs
    ctransext=avs
    netpbm=no

elif test -r /usr/local/apps/ncarg-4.2.2-pgi/lib/libncarg.a; then
#   ctransarg=avs
#   ctransext=x
    ctransarg=sun
    ctransext=sun
    netpbm=yes

elif test -r /opt/ncarg/bin/ctrans; then
    ctransarg=sun
    ctransext=sun
    netpbm=yes

elif test "$PLATFORM" = "AIX"; then
    ctransarg=avs
    ctransext=x

elif test "$PLATFORM" = "HP-UX"; then
    ctransarg=sgi
    ctransext=sgi

else
    ctransarg=sgi
    ctransext=sgi
#   ctransarg=sun
#   ctransext=sun

fi

echo "ctransarg=$ctransarg"
echo "ctransext=$ctransext"
echo "netpbm=$netpbm"

if test -r /opt/ncarg/bin/ctrans; then
    date
    echo "Running ctrans and netpbm programs to make gmeta_$prod.gif file on JET"
    COMFILE=$LAPS_DATA_ROOT/myctrans
    echo "COMFILE = $COMFILE"
    echo "#!/bin/csh" > $COMFILE                     
    echo "cd $LAPS_DATA_ROOT/lapsprd/www/anal2d; setenv NCARG_ROOT /opt/bin; ctrans -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta | rasttopnm | ppmtogif > $SCRATCH_DIR/gmeta_$prod.gif" >> $COMFILE                            
    chmod 775 $COMFILE                     
    cat $COMFILE                   
    echo "ssh v00 $COMFILE"                          
          ssh v00 $COMFILE                              
    date

#   For testing only...
    cp gmeta gmeta_$prod

elif test "$netpbm" = "yes"; then 
    if test -r /whomenull; then
        CTRANS=/opt/ncl/5.1.0_bin/bin/ctrans
#   elif test -r /usr/local/ncarg; then
#       CTRANS=/usr/local/ncarg/bin/ctrans
    elif test -r /usr/local/apps/ncarg-4.3.1.LINUX9; then
        CTRANS=/usr/local/apps/ncarg-4.3.1.LINUX9/bin/ctrans
    elif test -r /usr/local/ncarg-5.0.0-pgi-64-SLES10; then
        CTRANS=/usr/local/ncarg-5.0.0-pgi-64-SLES10/bin/ctrans
    elif test -r /usr/local/ncarg-5.0.0; then
        CTRANS=/usr/local/ncarg-5.0.0/bin/ctrans
    else
        CTRANS=$NCARG_ROOT/bin/ctrans
    fi

    date
    echo "Running $CTRANS and netpbm programs to make gmeta_$prod.gif file"
    which rasttopnm
    which ppmtogif
    echo "$CTRANS -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta | rasttopnm | ppmtogif > $SCRATCH_DIR/gmeta_$prod.gif"
          $CTRANS -verbose -d sun -window $WINDOW -resolution $RESOLUTION gmeta | rasttopnm | ppmtogif > $SCRATCH_DIR/gmeta_$prod.gif
    ls -l $SCRATCH_DIR/gmeta_$prod.gif
    date

else
    echo "Running $NCARG_ROOT/bin/ctrans to make gmeta_$prod.$ctransext file"

    $NCARG_ROOT/bin/ctrans -verbose -d $ctransarg -window $WINDOW -resolution $RESOLUTION gmeta > $SCRATCH_DIR/gmeta_$prod.$ctransext

    ls -l $SCRATCH_DIR/gmeta_$prod.$ctransext
    date
    echo "Running gmeta.$ctransext to GIF converter $LAPS_ETC/laps_prod.sh on $MACHINE"
    $LAPS_ETC/laps_prod.sh $prod $SCRATCH_DIR $SCRATCH_DIR $ctransext

    #Cleanup
    'rm' -f $SCRATCH_DIR/gmeta_$prod.$ctransext
    echo " "
    date
    echo " "

fi

chmod 666 $SCRATCH_DIR/gmeta_$prod.gif

#Copy and link the output GIF file
mkdir -p                                                               $WWW_DIR/anal2d/archive/$prod
mkdir -p                                                               $WWW_DIR/anal2d/recent/$prod
#mkdir -p                                                              $WWW_DIR/anal2d/loop/$prod

$CP $SCRATCH_DIR/gmeta_$prod.gif                                    $WWW_DIR/anal2d/archive/$prod/$datetime.gif

cd $WWW_DIR/anal2d/recent

rm -f                                                                gmeta_$prod$uscore$latest.gif
#ln -s ../archive/$prod/$datetime.gif                                gmeta_$prod$uscore$latest.gif
cp ../archive/$prod/$datetime.gif                                    gmeta_$prod$uscore$latest.gif

rm -f                                                                $prod/gmeta_$prod$uscore$latest.gif
#ln -s ../../archive/$prod/$datetime.gif                             $prod/gmeta_$prod$uscore$latest.gif
cp        ../archive/$prod/$datetime.gif                             $prod/gmeta_$prod$uscore$latest.gif

rm -f                                                                gmeta_$prod$uscore$utc_hhmm.gif
#ln -s ../archive/$prod/$datetime.gif                                gmeta_$prod$uscore$utc_hhmm.gif

rm -f                                                                $prod/gmeta_$prod$uscore$utc_hhmm.gif
#ln -s ../../archive/$prod/$datetime.gif                             $prod/gmeta_$prod$uscore$utc_hhmm.gif
cp        ../archive/$prod/$datetime.gif                             $prod/gmeta_$prod$uscore$utc_hhmm.gif

#ln -s $WWW_DIR/anal2d/archive/$prod/$datetime.gif                   $WWW_DIR/anal2d/loop/$prod/$datetime.gif

#Make the soft link for Web display of directories
cd    $WWW_DIR/anal2d/archive/$prod
rm -f                                                                  80.htaccess

#if test -r $WWW_DIR/../../utilities/public-current-directory.htaccess; then
#     ln -s $WWW_DIR/../../utilities/public-current-directory.htaccess       80.htaccess
#fi

#Create the files.txt script with the equivalent of files.sh
cd  $WWW_DIR/anal2d/loop/$prod
/bin/ls -1r *.gif > files.txt  

#Link to the files.cgi script
cd    $WWW_DIR/anal2d/loop/$prod

if test -r $WWW_DIR/../../looper/files.cgi; then
     rm -f files.cgi
     ln -s $WWW_DIR/../../looper/files.cgi files.cgi
else # this will use the supplied files.cgi
     rm -f files.cgi
     ln -s $LAPS_ETC/www/anal2d/files.cgi files.cgi
fi

echo "Finish laps_gifs_sub.sh..."
