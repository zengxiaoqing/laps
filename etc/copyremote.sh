#!/bin/sh

# Argument 1: filename

# Argument 2: Destination directory

# Argument 3: Copy method (bbcp, exchange, rsync)

# .................................................................................

FILENAME=$1
DESTDIR=$2
TEMPFILE=$FILENAME.$$
TEMPDIR=/exchange/tmp/fab

RSH=--rsh=ssh
RSYNCARGS=-rlptgvvz
REMOTE_NODE=pinky.fsl.noaa.gov

if test "$3" == "bbcp"; then
    echo "bbcp -s 8 -w 1M -P 5 -V $FILENAME $REMOTE_NODE:$DESTDIR"
          bbcp -s 8 -w 1M -P 5 -V $FILENAME $REMOTE_NODE:$DESTDIR

elif test "$3" == "exchange"; then
    echo "scp -r $FILENAME jetscp.rdhpcs.noaa.gov:$TEMPDIR/$TEMPFILE" 
          scp -r $FILENAME jetscp.rdhpcs.noaa.gov:$TEMPDIR/$TEMPFILE  

    if test -d /exchange/tmp/fab/$TEMPFILE; then # move individual files between the two directories
        echo " "
        echo "second hop of directory"
        echo "ssh $REMOTE_NODE mv $TEMPDIR/$TEMPFILE/* $DESTDIR/$FILENAME"
              ssh $REMOTE_NODE mv $TEMPDIR/$TEMPFILE/* $DESTDIR/$FILENAME 
    else
        echo "second hop of file"
        echo "ssh $REMOTE_NODE mv $TEMPDIR/$TEMPFILE $DESTDIR/$FILENAME"
              ssh $REMOTE_NODE mv $TEMPDIR/$TEMPFILE $DESTDIR/$FILENAME 
    fi

elif test "$3" == "rsync"; then
    echo "rsync $RSH $RSYNCARGS $FILENAME $REMOTE_NODE:$DESTDIR"                                    
          rsync $RSH $RSYNCARGS $FILENAME $REMOTE_NODE:$DESTDIR                                     

else
    echo "Usage: copyremote.sh [filename] [destination directory] [method]"

fi
