#!/bin/sh

file=$1
n_waits=$2
waittime=$3

n=0

while [ $n != $n_waits ]
    do
    n=`expr $n + 1`                                                     

    if [ -f $file ]
        then
        echo $file exists...
        n=$n_waits
    else
        echo $file does not exist $n...
        sleep $waittime
    fi

done

