#!/bin/csh -x

#This script sets up the pre-compiled tar file and localizes it at CWB
#Command line argument "p" does the configure/install for the pre-compiled tar file
#Command line argument "w" does the localization for the particular window

set arg1 = np
set arg2 = nw

foreach i ($argv)
    echo $i
    if ($i == p) then
        set arg1 = p
    endif

    if ($i == w) then
        set arg2 = w
    endif
end

#Primary LAPS environment variables
setenv LAPS_SRC_ROOT      /pj/fsl/albers/src/laps-0-8-19
setenv LAPSINSTALLROOT    /pj/fsl/albers/data_disk/laps
setenv LAPS_DATA_ROOT     /pj/fsl/albers/data_disk/laps/data
setenv TEMPLATEDIR        /pj/fsl/albers/src/template
setenv NEWPERL            /pj/fsl/albers/bin/perl

#This section does what "configure" would normally do for a precompiled tar file
if ($arg1 == p) then
    echo "Configuring precompiled scripts in $LAPSINSTALLROOT/etc"
    cd $LAPSINSTALLROOT/etc

#   Edit PERL paths

#   setenv OLDPERL      /usr/local/bin/perl
#   fix_net $OLDPERL    $NEWPERL   *.pl 
   
#   setenv OLDPERL      /usr/local/apps/perl5/bin/perl
#   fix_net $OLDPERL    $NEWPERL   *.pl 

#   setenv OLDPERL      /usr/local/perl5/bin/perl
#   fix_net $OLDPERL    $NEWPERL   *.pl 

#   setenv OLDPERL      /usr/local/perl/bin/perl
#   grep $OLDPERL       *.pl
#   fix_net $OLDPERL    $NEWPERL   *.pl 

#   setenv OLDPERL      /usr/nfs/bin/perl
#   grep $OLDPERL       *.pl
#   fix_net $OLDPERL    $NEWPERL   *.pl 

    setenv NEWNETCDF    /usr/local/netcdf/bin

#   Edit NetCDF paths
#   setenv OLDNETCDF    /usr/local/apps/netcdf/bin
#   fix_net $OLDNETCDF   $NEWNETCDF                 *.pl *.pm; echo " "

    foreach file (*.pl)
        if (-e $file.in) then
            cp $file.in $file
            fix_net @PERL@                   $NEWPERL            $file
            fix_net @NETCDF@                 $NEWNETCDF          $file
            fix_net @prefix@                 $LAPSINSTALLROOT    $file
            fix_net @datadir@                $LAPS_DATA_ROOT     $file
            fix_net @top_srcdir@             $LAPS_SRC_ROOT      $file
            fix_net @CSH@                    /bin/csh            $file

#           This edit is specific to Taiwan's @INC problem though it probably doesn't do much harm in other installations
            fix_net \$EXECUTABLE_NAME        $NEWPERL            $file

        endif
    end

else
    echo "Skipping configure/install step for precompiled tar file"

endif


#We assume that window_laps_rt.pl will know what to do if $LAPS_SRC_ROOT = $LAPS_DATA_ROOT
#                                                      just do localization with no templates
#
#We will see if window_laps_rt.pl will know what to do if $LAPS_SRC_ROOT != $LAPS_DATA_ROOT 
#                                                      and there is no $LAPS_SRC_ROOT/data/static/*.nl
#                                                      and there is no $LAPS_SRC_ROOT/data/cdl/*.cdl
#                                                      (e.g. precompiled tar file setup at CWB with split tree)
if ($arg2 == w) then
    echo "Starting the window setup step"

    if ($LAPS_SRC_ROOT/data != $LAPS_DATA_ROOT) then

        echo "We have a split directory tree"

        if (! -e $LAPS_SRC_ROOT/data) then
            echo "Warning: data directory not available in LAPS_SRC_ROOT, making it"
            mkdir $LAPS_SRC_ROOT/data
        endif

        if (! -e $LAPS_SRC_ROOT/data/static) then
            echo "Warning: static directory not available in LAPS_SRC_ROOT, making it"
            mkdir $LAPS_SRC_ROOT/data/static
        endif

        if (! -e $LAPS_SRC_ROOT/data/static/background.nl || ! -e $LAPS_SRC_ROOT/data/static/nest7grid.parms) then
            echo "Warning: namelists not available in LAPS_SRC_ROOT"

            if (! -e $LAPS_DATA_ROOT/data/static/background.nl || ! -e $LAPS_DATA_ROOT/data/static/nest7grid.parms) then
                echo "Error: namelists also not available in LAPS_DATA_ROOT, exit"
                exit
            endif

            echo "Doing a precopy"
            cp $LAPS_DATA_ROOT/static/*.nl    $LAPS_SRC_ROOT/data/static
            cp $LAPS_DATA_ROOT/static/*.parms $LAPS_SRC_ROOT/data/static
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

    $NEWPERL $LAPSINSTALLROOT/etc/window_laps_rt.pl -cf -t$TEMPLATEDIR -s$LAPS_SRC_ROOT -i$LAPSINSTALLROOT -d$LAPS_DATA_ROOT 

#   cd $LAPS_SRC_ROOT
#   make DATAROOT=$LAPS_DATA_ROOT mkdatadirs 
  
    echo " "
    echo "Returned from window_laps_rt.pl"

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

#   Now we can edit the cronfile
    echo " "
    echo "Setting up a cronfile in $LAPS_DATA_ROOT/cronfile..."
    cp $LAPSINSTALLROOT/util/cronfile.in $LAPS_DATA_ROOT/cronfile
    cd $LAPS_DATA_ROOT
    fix_net @PERL@    $NEWPERL         cronfile
    fix_net @prefix@  $LAPSINSTALLROOT cronfile
    fix_net @datadir@ $LAPS_DATA_ROOT  cronfile

else
    echo "Skipping window setup step"

endif


