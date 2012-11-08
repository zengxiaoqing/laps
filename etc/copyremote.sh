#!/bin/sh

# Argument 1: filename (individual file, subdirectory, or file list)

# Argument 2: Destination directory

# Argument 3: Copy method (bbcp, exchange, rsync, exchange_tar, exchange_tarfilelist)

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

elif test "$3" == "exchange_tar"; then
    if test -d $FILENAME; then                                                       
        echo " "
        echo "first hop of directory"
        cd $FILENAME
        echo 'tar cvf - . | ssh jetscp.rdhpcs.noaa.gov "mkdir -p $TEMPDIR/$TEMPFILE; cd $TEMPDIR/$TEMPFILE; tar xpvf -"'
              tar cvf - . | ssh jetscp.rdhpcs.noaa.gov "mkdir -p $TEMPDIR/$TEMPFILE; cd $TEMPDIR/$TEMPFILE; tar xpvf -"
    else
        echo " "
        echo "first hop of file"        
        echo "scp -r $FILENAME jetscp.rdhpcs.noaa.gov:$TEMPDIR/$TEMPFILE" 
              scp -r $FILENAME jetscp.rdhpcs.noaa.gov:$TEMPDIR/$TEMPFILE  
    fi

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

elif test "$3" == "exchange_tarfilelist"; then
    if test -e $FILENAME; then                                                       
        echo " "
        echo "make local tar file from filelist"
        echo  tar -T $FILENAME -cvf $TEMPFILE.tar    
              tar -T $FILENAME -cvf $TEMPFILE.tar                                                                         
        ls -l $TEMPFILE.tar

        echo  gzip $TEMPFILE.tar 
              gzip $TEMPFILE.tar 

        ls -l $TEMPFILE.tar.gz

        date -u

        echo " "
        echo "first hop of tar file"        
        echo  scp -r $TEMPFILE.tar.gz jetscp.rdhpcs.noaa.gov:$TEMPDIR              
              scp -r $TEMPFILE.tar.gz jetscp.rdhpcs.noaa.gov:$TEMPDIR               

        date -u

        TEMPFILE2=`basename $TEMPFILE`                        

        echo " "
        echo "zcat from $TEMPDIR into second directory"
        echo  ssh $REMOTE_NODE "cd $DESTDIR; zcat $TEMPDIR/$TEMPFILE2.tar.gz | tar -xvf -"
              ssh $REMOTE_NODE "cd $DESTDIR; zcat $TEMPDIR/$TEMPFILE2.tar.gz | tar -xvf -"
       
        date -u

    fi

    exit

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

elif test "$3" == "exchange_tar_filelist2"; then
    echo " "
    echo "first hop of file list"
    echo "tar -T $FILENAME -cvf - . | ssh jetscp.rdhpcs.noaa.gov mkdir -p $TEMPDIR/$TEMPFILE; cd $TEMPDIR/$TEMPFILE; tar xpvf -"
          tar -T $FILENAME -cvf - . | ssh jetscp.rdhpcs.noaa.gov mkdir -p $TEMPDIR/$TEMPFILE; cd $TEMPDIR/$TEMPFILE; tar xpvf -

elif test "$3" == "rsync"; then
    echo "rsync $RSH $RSYNCARGS $FILENAME $REMOTE_NODE:$DESTDIR"                                    
          rsync $RSH $RSYNCARGS $FILENAME $REMOTE_NODE:$DESTDIR                                     

else
    echo "Usage: copyremote.sh [filename] [destination directory] [method]"

fi
