#!/bin/csh 

echo " "
echo "Starting install_laps.csh...$6"

#This script sets up the pre-compiled tar file and localizes it at CWB
#Usage: install_laps.csh $LAPS_SRC_ROOT $LAPSINSTALLROOT $LAPS_DATA_ROOT $TEMPLATE `which perl` [a/w/p]

setenv LAPS_SRC_ROOT      $1 # Top level LAPS directory above 'src'
setenv LAPSINSTALLROOT    $2
setenv LAPS_DATA_ROOT     $3

#The TEMPLATEDIR is needed except for the 'p' option, in which case you can set it 
#to a dummy location such as '/dev/null'.
setenv TEMPLATEDIR        $4

#Note that this can normally be set to 'perl'. The full path should be used if
#the .cshrc does not have perl in its $path. The one other unusual time for
#passing in the full path be when a precompiled tar file is used and we are 
#using a special perl script to circumvent @INC problems (e.g. at CWB).
setenv NEWPERL            $5

# Argument $6
#    a - do entire install (equivalent to 'p' followed by 'w')
#    p - does binary (configure/install) setup only
#    w - does window localization only

#if ($i == p) then
#    set arg1 = p
#endif

#if ($i == w) then
#    set arg2 = w
#endif

#Primary LAPS environment variables
#setenv TEMPLATEDIR        /pj/fsl/albers/src/template
#setenv NEWPERL            /pj/fsl/albers/bin/perl
#setenv NEWPERL            `which perl`

if ($7 == intel) then
  echo "module switch pgi intel"
  module switch pgi intel
  which ifort
endif

#Configure/install?
if ($6 != w) then

# Precompiled tar file?
  if (! -e $LAPS_SRC_ROOT/Makefile) then

#   This section does what "configure" would normally do on a precompiled tar file,
#   (i.e. a substitution for section 2.2.2 in the README)

    echo "Makefile not present, assuming precompiled tar file"

    echo "Configuring precompiled scripts in $LAPSINSTALLROOT/etc and $LAPSINSTALLROOT/util"

    if(! -e $LAPS_SRC_ROOT/etc) then
        echo "ERROR: $LAPS_SRC_ROOT/etc is not present"
        exit
    endif

    cd $LAPS_SRC_ROOT

    setenv NEWNETCDF    /usr/local/netcdf

#   Edit paths

#   ls -l etc/*.pl util/cronfile

    foreach file (etc/*.pl util/cronfile)
        if (-e $file.in) then
            cp $file.in $file
            $LAPS_SRC_ROOT/util/fix_net @PERL@                   $NEWPERL            $file
            $LAPS_SRC_ROOT/util/fix_net @NETCDF@                 $NEWNETCDF          $file
            $LAPS_SRC_ROOT/util/fix_net @prefix@                 $LAPSINSTALLROOT    $file
#           $LAPS_SRC_ROOT/util/fix_net @datadir@                $LAPS_DATA_ROOT     $file
            $LAPS_SRC_ROOT/util/fix_net @top_srcdir@             $LAPS_SRC_ROOT      $file
            $LAPS_SRC_ROOT/util/fix_net @CSH@                    /bin/csh            $file

#           This edit is specific to Taiwan's @INC problem though it probably doesn't do much harm in other installations
            $LAPS_SRC_ROOT/util/fix_net \$EXECUTABLE_NAME        $NEWPERL            $file

        endif
    end

    ./util/fix_net .environs  .youll_never_find_this_file etc/fxa.pm

#   This section does what "make install" would normally do on a precompiled tar file,
#   (i.e. a substitution for section 2.2.4 in the README)

    mkdir -p $LAPSINSTALLROOT

    if ($LAPS_SRC_ROOT != $LAPSINSTALLROOT) then
        echo "We have a split directory tree with the INSTALLED binaries separated from SRC"

        cd $LAPS_SRC_ROOT

        if (-e bin) then
            echo "Moving $LAPS_SRC_ROOT/bin directory to $LAPSINSTALLROOT"
            mv bin  $LAPSINSTALLROOT
        endif

        echo "Copying $LAPS_SRC_ROOT/etc directory to $LAPSINSTALLROOT"
        cp -r etc  $LAPSINSTALLROOT

        echo "Copying $LAPS_SRC_ROOT/util directory to $LAPSINSTALLROOT"
        cp -r util  $LAPSINSTALLROOT
    endif

    chmod -R g+w $LAPSINSTALLROOT

  else
    echo "Makefile present, assuming a regular tar file"
    echo "We can proceed with configure/make if desired"

    cd $LAPS_SRC_ROOT

    echo " "
    echo "make"
    make >& make.out 
    ls -l make.out 

    echo " "
    echo "make install"
    make install >& make_install.out 
    ls -l make_install.out 

    echo " "
    echo "make install_lapsplot"
    make install_lapsplot >& make_install_lapsplot.out 
    ls -l make_install_lapsplot.out 

#   echo " "
#   echo "make install_wfopost"
#   make install_wfopost >& make_install_wfopost.out 
#   ls -l make_install_wfopost.out 

    echo " "
    echo "make debug"
    make debug >& make_debug.out 
    ls -l make_debug.out

    echo " "
    echo "chmod"
    chmod -R g+w $LAPS_SRC_ROOT

    cd $LAPSINSTALLROOT

    uname -p >& uname.out
    uname -n >> uname.out

  endif

else
  echo "Skipping configure/install step"

endif


#We assume that window_domain_rt.pl will know what to do if $LAPS_SRC_ROOT = $LAPS_DATA_ROOT
#                                                      just do localization with no templates
#
#We will see if window_domain_rt.pl will know what to do if $LAPS_SRC_ROOT != $LAPS_DATA_ROOT 
#                                                      and there is no $LAPS_SRC_ROOT/data/static/*.nl
#                                                      and there is no $LAPS_SRC_ROOT/data/cdl/*.cdl
#                                                      (e.g. precompiled tar file setup at CWB with split tree)
#See the README section 2.2.6.2
if ($6 != p) then
    echo "Starting the window setup step"

    if ($LAPS_SRC_ROOT/data != $LAPS_DATA_ROOT) then

        echo "We have a split directory tree with the DATA separated from SRC"

        if (! -e $LAPS_SRC_ROOT/data) then
            echo "Warning: data directory not available in LAPS_SRC_ROOT, making it"
            mkdir $LAPS_SRC_ROOT/data
        endif

        if (! -e $LAPS_SRC_ROOT/data/static) then
            echo "Warning: static directory not available in LAPS_SRC_ROOT, making it"
            mkdir $LAPS_SRC_ROOT/data/static
        endif

        if (! -e $LAPS_SRC_ROOT/data/static/background.nl || ! -e $LAPS_SRC_ROOT/data/static/nest7grid.parms) then
            echo "Error: namelists not available in LAPS_SRC_ROOT"
            exit

#           if (! -e $LAPS_DATA_ROOT/data/static/background.nl || ! -e $LAPS_DATA_ROOT/data/static/nest7grid.parms) then
#               echo "Error: namelists also not available in LAPS_DATA_ROOT, exit"
#               exit
#           endif

#           echo "Doing a precopy of namelists from $LAPS_DATA_ROOT to $LAPS_SRC_ROOT"
#           cp $LAPS_DATA_ROOT/static/*.nl    $LAPS_SRC_ROOT/data/static
#           cp $LAPS_DATA_ROOT/static/*.parms $LAPS_SRC_ROOT/data/static
        endif

        if (! -e $LAPS_SRC_ROOT/data/cdl/lga.cdl) then
            echo "Warning: cdls not available in LAPS_SRC_ROOT"

            if (! -e $LAPS_DATA_ROOT/cdl/lga.cdl) then
                echo "Error: cdls also not available in LAPS_DATA_ROOT, exit"
                exit
            endif

            echo "Doing a precopy"
            cp $LAPS_DATA_ROOT/cdl/*.cdl $LAPS_SRC_ROOT/data/cdl
        endif
    endif

    rm -rf $LAPS_DATA_ROOT/*_save

    if (-e /home/oplapb/jet) then
        setenv QSPN "-q usfsfire"
    else if (-e /home/oplapb/ijet) then
        setenv QSPN "-q usfsfire"
    else
        setenv QSPN " "
    endif

    if (-e $LAPS_DATA_ROOT/lapsprd) then
        setenv config_domain f
        echo " "
        echo "Calling window_domain_rt.pl, config_domain = f"
        $NEWPERL $LAPSINSTALLROOT/etc/window_domain_rt.pl    -t $TEMPLATEDIR -s $LAPS_SRC_ROOT -i $LAPSINSTALLROOT -d $LAPS_DATA_ROOT -w laps $QSPN
    else
        setenv config_domain t
        echo " "
        echo "Calling window_domain_rt.pl, config_domain = t"
        $NEWPERL $LAPSINSTALLROOT/etc/window_domain_rt.pl -c -t $TEMPLATEDIR -s $LAPS_SRC_ROOT -i $LAPSINSTALLROOT -d $LAPS_DATA_ROOT -w laps $QSPN
    endif

#   cd $LAPS_SRC_ROOT
#   make DATAROOT=$LAPS_DATA_ROOT mkdatadirs 
  
    echo " "
    echo "Returned from window_domain_rt.pl"

    echo " "
    echo "checking LAPS_SRC_ROOT subdirs (static & cdl)..."
    ls -l $LAPS_SRC_ROOT/data/static
    ls -l $LAPS_SRC_ROOT/data/cdl

    echo " "
    echo "checking LAPS_DATA_ROOT lapsprd subdirs..."
    ls -l $LAPS_DATA_ROOT/lapsprd

    echo " "
    echo "checking static.nest7grid..."
    ls -l $LAPS_DATA_ROOT/static/static.nest7grid

#   Now we can edit the cronfile (README section 2.4)
    echo " "
    echo "Setting up a cronfile in $LAPS_DATA_ROOT..."
    cd $LAPSINSTALLROOT/etc
    $NEWPERL cronfile.pl --installroot=$LAPSINSTALLROOT --dataroot=$LAPS_DATA_ROOT

else
    echo "Skipping window setup step"

endif

if (-e /home/oplapb/lapsplot.exe) then
    echo "Copy pregenerated version of lapsplot.exe"
    cp /home/oplapb/lapsplot.exe $LAPSINSTALLROOT/bin
endif

echo "End of install_laps.csh..."


