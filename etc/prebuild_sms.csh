#!/bin/csh

setenv $LAPS_SRC_ROOT $1
setenv srcdir $2

cd $LAPS_SRC_ROOT/$srcdir
pwd

setenv SMS /home/schaffer/sms.steve

cp $LAPS_SRC_ROOT/src/include/directives.inc .

foreach filename (*.f)
    setenv file `echo $filename | cut -f 1 -d .`
    ls $file.f

    $SMS/bin/ppp --comment --V=1 --header directives.inc
    cat $file | gcc -E -P -traditional - > $file_cpp.f
    $SMS/bin/ppp --includepath=$SMS/include --Fcommon=directives.inc --comment --V=1 $file_cpp.f
    mv $file_cpp_sms.f $file.f

end


