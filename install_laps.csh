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
#setenv NEWPERL            /pj/fsl/albers/bin/perl
setenv NEWPERL            `which perl`

#Configure/install a precompiled tar file
if ($arg1 == p) then

#   This section does what "configure" would normally do on a precompiled tar file,
#   (i.e. a substitution for section 2.2.2 in the README)

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
            $LAPS_SRC_ROOT/util/fix_net @PERL@                   $NEWPERL            $file
            $LAPS_SRC_ROOT/util/fix_net @NETCDF@                 $NEWNETCDF          $file
            $LAPS_SRC_ROOT/util/fix_net @prefix@                 $LAPSINSTALLROOT    $file
            $LAPS_SRC_ROOT/util/fix_net @datadir@                $LAPS_DATA_ROOT     $file
            $LAPS_SRC_ROOT/util/fix_net @top_srcdir@             $LAPS_SRC_ROOT      $file
            $LAPS_SRC_ROOT/util/fix_net @CSH@                    /bin/csh            $file

#           This edit is specific to Taiwan's @INC problem though it probably doesn't do much harm in other installations
            $LAPS_SRC_ROOT/util/fix_net \$EXECUTABLE_NAME        $NEWPERL            $file

        endif
    end

#   This section does what "make install" would normally do on a precompiled tar file,
#   (i.e. a substitution for section 2.2.4 in the README)

    mkdir -p $LAPSINSTALLROOT

    if ($LAPS_SRC_ROOT != $LAPSINSTALLROOT) then
        echo "We have a split directory tree"

        cd $LAPS_SRC_ROOT

        if (-e bin) then
            echo "Moving $LAPS_SRC_ROOT/bin directory to $LAPSINSTALLROOT"
            mv bin  $LAPSINSTALLROOT
        endif

        if (-e etc) then
            echo "Moving $LAPS_SRC_ROOT/etc directory to $LAPSINSTALLROOT"
            mv etc  $LAPSINSTALLROOT
        endif

        if (-e util) then
            echo "Moving $LAPS_SRC_ROOT/util directory to $LAPSINSTALLROOT"
            mv util  $LAPSINSTALLROOT
        endif
    endif

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
#See the README section 2.2.6.2
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

#   Now we can edit the cronfile (README section 2.4)
    echo " "
    echo "Setting up a cronfile in $LAPS_DATA_ROOT/cronfile..."
    cp $LAPSINSTALLROOT/util/cronfile.in $LAPS_DATA_ROOT/cronfile
    cd $LAPS_DATA_ROOT
    $LAPS_SRC_ROOT/util/fix_net @PERL@    $NEWPERL         cronfile
    $LAPS_SRC_ROOT/util/fix_net @prefix@  $LAPSINSTALLROOT cronfile
    $LAPS_SRC_ROOT/util/fix_net @datadir@ $LAPS_DATA_ROOT  cronfile

else
    echo "Skipping window setup step"

endif


