#!/bin/csh

setenv $LAPS_SRC_ROOT $1
setenv srcdir $2
setenv SMS $3

echo "begin prebuild_sms.csh - SMS environment variable is $SMS"

cd $LAPS_SRC_ROOT/$srcdir
pwd


cp $LAPS_SRC_ROOT/src/include/directives.inc .

foreach filename (*.f)
    setenv file `echo $filename | cut -f 1 -d .`
    setenv nsms `grep -i csms $file.f | wc -l`
    echo "$file.f $nsms"

#   Translate with PPP only those files with SMS directives within
    if ($nsms != 0) then
        $SMS/bin/ppp --debug --comment --V=1 --header directives.inc
        cat $file.f | gcc -E -P -traditional - > $file\_cpp.f
        $SMS/bin/ppp --debug --includepath=$SMS/include --includepath=$LAPS_SRC_ROOT/src/include --Fcommon=directives.inc --comment --V=1 $file\_cpp.f
        mv $file\_cpp\_sms.f $file.f
    endif

end


