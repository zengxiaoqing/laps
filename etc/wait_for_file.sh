#!/bin/sh

n=1

file=$1

while [ $n != $2 ]
    do
    n=`expr $n + 1`                                                     

    if [ -f $file ]
        then
        echo $file exists...
        sleep 60
    else
        echo $file does not exist...
        n=$2
    fi

done

