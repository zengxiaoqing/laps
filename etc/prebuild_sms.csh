#!/bin/csh

setenv $LAPS_SRC_ROOT $1
setenv srcdir $2

echo "begin prebuild_sms.csh..."

cd $LAPS_SRC_ROOT/$srcdir
pwd

setenv SMS /home/schaffer/sms.steve

cp $LAPS_SRC_ROOT/src/include/directives.inc .

foreach filename (*.f)
    setenv file `echo $filename | cut -f 1 -d .`
    ls $file.f

    $SMS/bin/ppp --comment --V=1 --header directives.inc
    cat $file.f | gcc -E -P -traditional - > $file\_cpp.f
    $SMS/bin/ppp --includepath=$SMS/include --Fcommon=directives.inc --comment --V=1 $file\_cpp.f
    mv $file\_cpp\_sms.f $file.f

end


