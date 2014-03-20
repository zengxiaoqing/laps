#!/bin/sh

# Argument 1: filename (individual file, subdirectory, or file list)

# Argument 2: Destination directory (without remote node)

# Argument 3: Copy method (bbcp, exchange, rsync, exchange_tar, exchange_tarfilelist, 
#                          exchange_tarfilelist_z, exchange_filelist, exchange_zip, exchange_zip2)

# .................................................................................

FILENAME=$1
DESTDIR=$2
TEMPFILE=$FILENAME.$$
TEMPDIR=/exchange/tmp/fab/zeus

RSH=--rsh=ssh
RSYNCARGS=-rlptgvvz
RSYNCARGS_NOCOMPRESS=-rlptgvv 
REMOTE_NODE=clank.fsl.noaa.gov
REMOTE_NODE_MV=dlaps-ms1.fsl.noaa.gov

echo "copy_remote.sh..."

if test "$3" == "bbcp"; then
    echo "bbcp -s 32 -w 1M -P 5 -V $FILENAME $REMOTE_NODE:$DESTDIR"
          bbcp -s 32 -w 1M -P 5 -V $FILENAME $REMOTE_NODE:$DESTDIR

elif test "$3" == "exchange"; then # best for subdirectories
    echo "scp -o ConnectTimeout=200 -r $FILENAME jetscp.rdhpcs.noaa.gov:$TEMPDIR/$TEMPFILE" 
          scp -o ConnectTimeout=200 -r $FILENAME jetscp.rdhpcs.noaa.gov:$TEMPDIR/$TEMPFILE  

    if test -d $FILENAME; then # move individual files between the two directories
        echo " "
        date -u
        echo "second hop of directory"
        echo "ssh $REMOTE_NODE_MV mv $TEMPDIR/$TEMPFILE/* $DESTDIR/$FILENAME"
              ssh $REMOTE_NODE_MV mv $TEMPDIR/$TEMPFILE/* $DESTDIR/$FILENAME 
    else
        date -u
        echo "second hop of file"
        echo "ssh $REMOTE_NODE_MV mv $TEMPDIR/$TEMPFILE $DESTDIR/$FILENAME"
              ssh $REMOTE_NODE_MV mv $TEMPDIR/$TEMPFILE $DESTDIR/$FILENAME 
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
        echo "scp -o ConnectTimeout=200 -r $FILENAME jetscp.rdhpcs.noaa.gov:$TEMPDIR/$TEMPFILE" 
              scp -o ConnectTimeout=200 -r $FILENAME jetscp.rdhpcs.noaa.gov:$TEMPDIR/$TEMPFILE  
    fi

    if test -d $FILENAME; then # move individual files between the two directories
        echo " "
        date -u
        echo "second hop of directory"
        echo "ssh $REMOTE_NODE mv $TEMPDIR/$TEMPFILE/* $DESTDIR/$FILENAME"
              ssh $REMOTE_NODE mv $TEMPDIR/$TEMPFILE/* $DESTDIR/$FILENAME 
    else
        date -u
        echo "second hop of file"
        echo "ssh $REMOTE_NODE mv $TEMPDIR/$TEMPFILE $DESTDIR/$FILENAME"
              ssh $REMOTE_NODE mv $TEMPDIR/$TEMPFILE $DESTDIR/$FILENAME 
    fi

elif test "$3" == "exchange_filelist"; then
    BASEFILE=`basename $FILENAME` 
    TEMPPATH=`dirname  $FILENAME` 

    echo " "
    echo "first hop of filelist"        
    echo "ssh $REMOTE_NODE mkdir -p $TEMPDIR/$BASEFILE.$$"
          ssh $REMOTE_NODE mkdir -p $TEMPDIR/$BASEFILE.$$

#   echo "scp -o ConnectTimeout=200 `cat $FILENAME` jetscp.rdhpcs.noaa.gov:$TEMPDIR/$TEMPFILE" 
#         scp -o ConnectTimeout=200 `cat $FILENAME` jetscp.rdhpcs.noaa.gov:$TEMPDIR/$TEMPFILE

#   echo "tar -T $FILENAME cvf - | ssh jetscp.rdhpcs.noaa.gov tar xpvf - -C $TEMPDIR/$TEMPFILE"
#         tar -T $FILENAME cvf - | ssh jetscp.rdhpcs.noaa.gov tar xpvf - -C $TEMPDIR/$TEMPFILE

    cd $TEMPPATH/..

    echo "Tarring file list from the following directory: $cwd"
    pwd

    echo "tar -T $BASEFILE -cvf - | (ssh jetscp.rdhpcs.noaa.gov; cd $TEMPDIR/$BASEFILE.$$;  tar xfBp -)"
          tar -T $BASEFILE -cvf - | (ssh jetscp.rdhpcs.noaa.gov; cd $TEMPDIR/$BASEFILE.$$;  tar xfBp -)

    echo " "
    date -u
    echo "second hop of filelist"
    echo "ssh $REMOTE_NODE mv $TEMPDIR/$BASEFILE.$$/* $DESTDIR"
          ssh $REMOTE_NODE mv $TEMPDIR/$BASEFILE.$$/* $DESTDIR 

    echo "ssh $REMOTE_NODE rmdir $TEMPDIR/$TEMPFILE"
          ssh $REMOTE_NODE rmdir $TEMPDIR/$TEMPFILE 

elif test "$3" == "exchange_tarfilelist"; then
    if test -e $FILENAME; then                                                       
        echo " "
        echo "make local tar file from filelist"
        echo  tar -T $FILENAME -cvf $TEMPFILE.tar
              tar -T $FILENAME -cvf $TEMPFILE.tar
        ls -l $TEMPFILE.tar

        date -u

        echo " "
        echo "first hop of tar file"        
        echo  scp -o ConnectTimeout=200 -r $TEMPFILE.tar jetscp.rdhpcs.noaa.gov:$TEMPDIR              
              scp -o ConnectTimeout=200 -r $TEMPFILE.tar jetscp.rdhpcs.noaa.gov:$TEMPDIR               

        date -u

        TEMPFILE2=`basename $TEMPFILE`                        

        echo " "
        echo "untar from $TEMPDIR into second directory"
        echo  ssh $REMOTE_NODE "cd $DESTDIR; tar -xvf $TEMPDIR/$TEMPFILE2.tar; rm -f $TEMPDIR/$TEMPFILE2.tar"
              ssh $REMOTE_NODE "cd $DESTDIR; tar -xvf $TEMPDIR/$TEMPFILE2.tar; rm -f $TEMPDIR/$TEMPFILE2.tar"
       
        date -u

    fi

elif test "$3" == "exchange_tarfilelist_z"; then
    if test -e $FILENAME; then                                                       
        echo " "
        echo "make local tar file from filelist"
        echo  tar -T $FILENAME --gzip -cvf $TEMPFILE.tar.gz    
              tar -T $FILENAME --gzip -cvf $TEMPFILE.tar.gz                                                                         
        ls -l $TEMPFILE.tar.gz

        date -u

        echo " "
        echo "first hop of tar file"        
        echo  scp -o ConnectTimeout=200 -r $TEMPFILE.tar.gz jetscp.rdhpcs.noaa.gov:$TEMPDIR              
              scp -o ConnectTimeout=200 -r $TEMPFILE.tar.gz jetscp.rdhpcs.noaa.gov:$TEMPDIR               

        date -u

        TEMPFILE2=`basename $TEMPFILE`                        

        echo " "
        echo "zcat from $TEMPDIR into second directory"
        echo  ssh $REMOTE_NODE "cd $DESTDIR; zcat $TEMPDIR/$TEMPFILE2.tar.gz | tar -xvf -; rm -f $TEMPDIR/$TEMPFILE2.tar.gz"
              ssh $REMOTE_NODE "cd $DESTDIR; zcat $TEMPDIR/$TEMPFILE2.tar.gz | tar -xvf -; rm -f $TEMPDIR/$TEMPFILE2.tar.gz"
       
        date -u

    fi

elif test "$3" == "exchange_zip"; then
    if test -e $FILENAME; then                                                       
        echo $DESTDIR > destdirname

        if test ! -f destdirname; then
            echo "ERROR: unable to create destdirname"
            exit
        fi

        echo " "
        echo "make local zip file from filelist"
        echo  zip $TEMPFILE.z destdirname `cat $FILENAME`   
              zip $TEMPFILE.z destdirname `cat $FILENAME`                                                                         
        ls -l $TEMPFILE.z

        date -u

        echo " "
        echo "first hop of zip file"        
        echo  rsync $RSH $RSYNCARGS_NOCOMPRESS $TEMPFILE.z jetscp.rdhpcs.noaa.gov:$TEMPDIR              
              rsync $RSH $RSYNCARGS_NOCOMPRESS $TEMPFILE.z jetscp.rdhpcs.noaa.gov:$TEMPDIR               

        date -u

    fi

elif test "$3" == "exchange_zip2"; then # for files or file list
    if test -e $FILENAME; then                                                       
        echo $DESTDIR > destdirname

        if test ! -f destdirname; then
            echo "ERROR: unable to create destdirname"
            exit
        fi

        echo " "
        echo "make local zip file from filelist"
        echo  zip $TEMPFILE.z2 destdirname `cat $FILENAME`   
              zip $TEMPFILE.z2 destdirname `cat $FILENAME`                                                                         
        ls -l $TEMPFILE.z2

        date -u

        echo " "
        echo "first hop of zip file"        
        echo  scp -o ConnectTimeout=200 $TEMPFILE.z2 jetscp.rdhpcs.noaa.gov:$TEMPDIR              
              scp -o ConnectTimeout=200 $TEMPFILE.z2 jetscp.rdhpcs.noaa.gov:$TEMPDIR               

        echo "exit status (first scp try) is $?"

        if test "$?" != "0"; then 
            echo "error scp failed, trying again..." 
            echo  scp -o ConnectTimeout=200 $TEMPFILE.z2 jetscp.rdhpcs.noaa.gov:$TEMPDIR              
                  scp -o ConnectTimeout=200 $TEMPFILE.z2 jetscp.rdhpcs.noaa.gov:$TEMPDIR               
        fi

        echo "exit status (second scp try) is $?"

        date -u

        TEMPFILE2=`basename $TEMPFILE`                        

        echo " "
        echo "unzip from $TEMPDIR into second directory"
        echo  ssh $REMOTE_NODE "cd $DESTDIR; unzip -o $TEMPDIR/$TEMPFILE2.z2; rm -f $TEMPDIR/$TEMPFILE2.z2"
              ssh $REMOTE_NODE "cd $DESTDIR; unzip -o $TEMPDIR/$TEMPFILE2.z2; rm -f $TEMPDIR/$TEMPFILE2.z2"

        echo "exit status is $?"
       
        date -u

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
