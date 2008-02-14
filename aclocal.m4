dnl
dnl Autoconfig functions adapted to the NOAA FSL LAPS 
dnl based on version 2.12 of gnu autoconfig    Jim Edwards 1997
dnl
dnl AC_INIT_LAPS_NOTICE()
AC_DEFUN(AC_INIT_LAPS_NOTICE,
[# Guess values for system-dependent variables and create Makefiles.
# Generated automatically using autoconf version] AC_ACVERSION [
# Copyright (C) 1992, 93, 94, 95, 96 Free Software Foundation, Inc.
#
# This configure script is free software; the Free Software Foundation
# gives unlimited permission to copy, distribute and modify it.
#
# Modified for LAPS (Local Analysis and Prediction System) 
# By Jim Edwards                      NOAA Forecast Systems Laboratory 
# 1997
# 

# Defaults:
ac_help=
# ac_default_prefix=/usr/local/laps
# ac_default_prefix=`pwd`
AC_PREFIX_DEFAULT(`pwd`)
[#] Any additions from configure.in:])

dnl AC_INIT_LAPS_PARSE_ARGS()
AC_DEFUN(AC_INIT_LAPS_PARSE_ARGS,
[
# Initialize some variables set by options.
# The variables have the same names as the options, with
# dashes changed to underlines.
build=NONE
cache_file=./config.cache
exec_prefix=NONE
host=NONE
no_create=
nonopt=NONE
no_recursion=
prefix=NONE
datadir=NONE
program_prefix=NONE
program_suffix=NONE
program_transform_name=s,x,x,
silent=
site=
srcdir=
target=NONE
verbose=
x_includes=NONE
x_libraries=NONE
dnl Installation directory options.
dnl These are left unexpanded so users can "make install exec_prefix=/foo"
dnl and all the variables that are supposed to be based on exec_prefix
dnl by default will actually change.
dnl Use braces instead of parens because sh, perl, etc. also accept them.
bindir='${exec_prefix}/bin'
sbindir='${exec_prefix}/sbin'
libexecdir='${exec_prefix}/libexec'
datadir='${prefix}/share'
sysconfdir='${prefix}/etc'
sharedstatedir='${prefix}/com'
localstatedir='${prefix}/var'
libdir='${exec_prefix}/lib'
includedir='${prefix}/include'
oldincludedir='/usr/include'
infodir='${prefix}/info'
mandir='${prefix}/man'

# Initialize some other variables.
subdirs=
MFLAGS= MAKEFLAGS=
# Maximum number of lines to put in a shell here document.
ac_max_here_lines=12

ac_prev=
for ac_option
do

  # If the previous option needs an argument, assign it.
  if test -n "$ac_prev"; then
    eval "$ac_prev=\$ac_option"
    ac_prev=
    continue
  fi

  case "$ac_option" in
changequote(, )dnl
  -*=*) ac_optarg=`echo "$ac_option" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
changequote([, ])dnl
  *) ac_optarg= ;;
  esac

  # Accept the important Cygnus configure options, so we can diagnose typos.

  case "$ac_option" in

  -bindir | --bindir | --bindi | --bind | --bin | --bi)
    ac_prev=bindir ;;
  -bindir=* | --bindir=* | --bindi=* | --bind=* | --bin=* | --bi=*)
    bindir="$ac_optarg" ;;

  -build | --build | --buil | --bui | --bu)
    ac_prev=build ;;
  -build=* | --build=* | --buil=* | --bui=* | --bu=*)
    build="$ac_optarg" ;;

  -cache-file | --cache-file | --cache-fil | --cache-fi \
  | --cache-f | --cache- | --cache | --cach | --cac | --ca | --c)
    ac_prev=cache_file ;;
  -cache-file=* | --cache-file=* | --cache-fil=* | --cache-fi=* \
  | --cache-f=* | --cache-=* | --cache=* | --cach=* | --cac=* | --ca=* | --c=*)
    cache_file="$ac_optarg" ;;

  -datadir | --datadir | --datadi | --datad | --data | --dat | --da)
    ac_prev=datadir ;;
  -datadir=* | --datadir=* | --datadi=* | --datad=* | --data=* | --dat=* \
  | --da=*)
    datadir="$ac_optarg" ;;

  -disable-* | --disable-*)
    ac_feature=`echo $ac_option|sed -e 's/-*disable-//'`
    # Reject names that are not valid shell variable names.
changequote(, )dnl
    if test -n "`echo $ac_feature| sed 's/[-a-zA-Z0-9_]//g'`"; then
changequote([, ])dnl
      AC_MSG_ERROR($ac_feature: invalid feature name)
    fi
    ac_feature=`echo $ac_feature| sed 's/-/_/g'`
    eval "enable_${ac_feature}=no" ;;

  -enable-* | --enable-*)
    ac_feature=`echo $ac_option|sed -e 's/-*enable-//' -e 's/=.*//'`
    # Reject names that are not valid shell variable names.
changequote(, )dnl
    if test -n "`echo $ac_feature| sed 's/[-_a-zA-Z0-9]//g'`"; then
changequote([, ])dnl
      AC_MSG_ERROR($ac_feature: invalid feature name)
    fi
    ac_feature=`echo $ac_feature| sed 's/-/_/g'`
    case "$ac_option" in
      *=*) ;;
      *) ac_optarg=yes ;;
    esac
    eval "enable_${ac_feature}='$ac_optarg'" ;;

  -exec-prefix | --exec_prefix | --exec-prefix | --exec-prefi \
  | --exec-pref | --exec-pre | --exec-pr | --exec-p | --exec- \
  | --exec | --exe | --ex)
    ac_prev=exec_prefix ;;
  -exec-prefix=* | --exec_prefix=* | --exec-prefix=* | --exec-prefi=* \
  | --exec-pref=* | --exec-pre=* | --exec-pr=* | --exec-p=* | --exec-=* \
  | --exec=* | --exe=* | --ex=*)
    exec_prefix="$ac_optarg" ;;

  -gas | --gas | --ga | --g)
    # Obsolete; use --with-gas.
    with_gas=yes ;;

  -help | --help | --hel | --he)
    # Omit some internal or obsolete options to make the list less imposing.
    # This message is too long to be a string in the A/UX 3.1 sh.
    cat << EOF
changequote(, )dnl
Usage: configure [options] [host]
Options: [defaults in brackets after descriptions]
Configuration:
  --cache-file=FILE       cache test results in FILE
  --help                  print this message
  --no-create             do not create output files
  --quiet, --silent       do not print \`checking...' messages
  --version               print the version of autoconf that created configure
Compilers:
  --cc=CC                 name of the C compiler to use [guessed]
  --cflags=CFLAGS         options to the C compiler [preset based on arch]
  --fc=f90                name of the Fortran compiler to use [guessed]
  --fflags=FFLAGS         options to the Fortran compiler [preset based on arch]
  --optimization=OPT      optimization options passed to both C and fortran compilers  
  
Directory and file names:
  --prefix=PREFIX         install LAPS files in PREFIX
                          [$ac_default_prefix]
  --datadir=DIR           LAPS output data in DIR
                          [PREFIX/data]
EOF
    cat << EOF
Host type:
  --arch=ARCH             configure for building on ARCH [guessed]
Features and packages:
  --netcdf=NETCDF         path to netcdf package (required) [guessed]
changequote([, ])dnl
EOF
    if test -n "$ac_help"; then
      echo "--enable and --with options recognized:$ac_help"
    fi
    exit 0 ;;

  -arch | --arch | --arc | --ar | --a)
    ac_prev=arch ;;
  -arch=* | --arch=* | --arc=* | --ar=* | --a=*)
    arch="$ac_optarg" ;;
  
  -cc | --cc)
    ac_prec=cc ;;
  -cc=* | --cc=*)
    CC="$ac_optarg" ;;

  -cflags | --cflags | --cflag | --cfla | --cfl | --cf)
    ac_prec=cflags ;;
  -cflags=* | --cflags=* | --cflag=* | --cfla=* | --cfl=* | --cf=*)
    user_cflags="$ac_optarg" ;;

  -fc | --fc)
    ac_prec=fc ;;
  -fc=* | --fc=*)
    FC="$ac_optarg" ;;

  -fflags | --fflags | --fflag | --ffla | --ffl | --ff)
    ac_prec=fflags ;;
  -fflags=* | --fflags=* | --fflag=* | --ffla=* | --ffl=* | --ff=*)
    user_fflags="$ac_optarg" ;;

  -includedir | --includedir | --includedi | --included | --include \
  | --includ | --inclu | --incl | --inc)
    ac_prev=includedir ;;
  -includedir=* | --includedir=* | --includedi=* | --included=* | --include=* \
  | --includ=* | --inclu=* | --incl=* | --inc=*)
    includedir="$ac_optarg" ;;

  -infodir | --infodir | --infodi | --infod | --info | --inf)
    ac_prev=infodir ;;
  -infodir=* | --infodir=* | --infodi=* | --infod=* | --info=* | --inf=*)
    infodir="$ac_optarg" ;;

  -libdir | --libdir | --libdi | --libd)
    ac_prev=libdir ;;
  -libdir=* | --libdir=* | --libdi=* | --libd=*)
    libdir="$ac_optarg" ;;

  -libexecdir | --libexecdir | --libexecdi | --libexecd | --libexec \
  | --libexe | --libex | --libe)
    ac_prev=libexecdir ;;
  -libexecdir=* | --libexecdir=* | --libexecdi=* | --libexecd=* | --libexec=* \
  | --libexe=* | --libex=* | --libe=*)
    libexecdir="$ac_optarg" ;;

  -localstatedir | --localstatedir | --localstatedi | --localstated \
  | --localstate | --localstat | --localsta | --localst \
  | --locals | --local | --loca | --loc | --lo)
    ac_prev=localstatedir ;;
  -localstatedir=* | --localstatedir=* | --localstatedi=* | --localstated=* \
  | --localstate=* | --localstat=* | --localsta=* | --localst=* \
  | --locals=* | --local=* | --loca=* | --loc=* | --lo=*)
    localstatedir="$ac_optarg" ;;

  -mandir | --mandir | --mandi | --mand | --man | --ma | --m)
    ac_prev=mandir ;;
  -mandir=* | --mandir=* | --mandi=* | --mand=* | --man=* | --ma=* | --m=*)
    mandir="$ac_optarg" ;;

  -ncargraphics | --ncargraphics | --ncargraphic | --ncargraphi \
  | --ncargraph | --ncargrap | --ncargra | --ncargr | --ncarg \
  | --ncar | --nca | --nc)
    ac_prev=ncargraphics ;;
  -ncargraphics=* | --ncargraphics=* | --ncargraphic=* | --ncargraphi=* \
  | --ncargraph=* | --ncargrap=* | --ncargra=* | --ncargr=* | --ncarg=* \
  | --ncar=* | --nca=* | --nc=*)
    ncargraphics="$ac_optarg" ;;

  -netcdf | --netcdf | --netcd | --netc | --net | --ne)
    ac_prev=netcdf ;;
  -netcdf=* | --netcdf=* | --netcd=* | --netc=* | --net=* | --ne=*)
    netcdf="$ac_optarg" ;;

  -no-create | --no-create | --no-creat | --no-crea | --no-cre \
  | --no-cr | --no-c)
    no_create=yes ;;

  -no-recursion | --no-recursion | --no-recursio | --no-recursi \
  | --no-recurs | --no-recur | --no-recu | --no-rec | --no-re | --no-r)
    no_recursion=yes ;;

  -optimization | --optimization | --optimizatio | --optimizati \
  | --optimizat | --optimiza | --optimiz | --optimi | --optim \
  | --opti | --opt | --op | --o)
    ac_prev=optimization ;;
  -optimization=* | --optimization=* | --optimizatio=* | --optimizati=* \
  | --optimizat=* | --optimiza=* | --optimiz=* | --optimi=* | --optim=* \
  | --opti=* | --opt=* | --op=* | --o=*)
    user_optimize="$ac_optarg" ;;	

  -prefix | --prefix | --prefi | --pref | --pre | --pr | --p)
    ac_prev=prefix ;;
  -prefix=* | --prefix=* | --prefi=* | --pref=* | --pre=* | --pr=* | --p=*)
    prefix="$ac_optarg" ;;

  -program-prefix | --program-prefix | --program-prefi | --program-pref \
  | --program-pre | --program-pr | --program-p)
    ac_prev=program_prefix ;;
  -program-prefix=* | --program-prefix=* | --program-prefi=* \
  | --program-pref=* | --program-pre=* | --program-pr=* | --program-p=*)
    program_prefix="$ac_optarg" ;;

  -program-suffix | --program-suffix | --program-suffi | --program-suff \
  | --program-suf | --program-su | --program-s)
    ac_prev=program_suffix ;;
  -program-suffix=* | --program-suffix=* | --program-suffi=* \
  | --program-suff=* | --program-suf=* | --program-su=* | --program-s=*)
    program_suffix="$ac_optarg" ;;

  -program-transform-name | --program-transform-name \
  | --program-transform-nam | --program-transform-na \
  | --program-transform-n | --program-transform- \
  | --program-transform | --program-transfor \
  | --program-transfo | --program-transf \
  | --program-trans | --program-tran \
  | --progr-tra | --program-tr | --program-t)
    ac_prev=program_transform_name ;;
  -program-transform-name=* | --program-transform-name=* \
  | --program-transform-nam=* | --program-transform-na=* \
  | --program-transform-n=* | --program-transform-=* \
  | --program-transform=* | --program-transfor=* \
  | --program-transfo=* | --program-transf=* \
  | --program-trans=* | --program-tran=* \
  | --progr-tra=* | --program-tr=* | --program-t=*)
    program_transform_name="$ac_optarg" ;;

  -q | -quiet | --quiet | --quie | --qui | --qu | --q \
  | -silent | --silent | --silen | --sile | --sil)
    silent=yes ;;

  -sbindir | --sbindir | --sbindi | --sbind | --sbin | --sbi | --sb)
    ac_prev=sbindir ;;
  -sbindir=* | --sbindir=* | --sbindi=* | --sbind=* | --sbin=* \
  | --sbi=* | --sb=*)
    sbindir="$ac_optarg" ;;

  -sharedstatedir | --sharedstatedir | --sharedstatedi \
  | --sharedstated | --sharedstate | --sharedstat | --sharedsta \
  | --sharedst | --shareds | --shared | --share | --shar \
  | --sha | --sh)
    ac_prev=sharedstatedir ;;
  -sharedstatedir=* | --sharedstatedir=* | --sharedstatedi=* \
  | --sharedstated=* | --sharedstate=* | --sharedstat=* | --sharedsta=* \
  | --sharedst=* | --shareds=* | --shared=* | --share=* | --shar=* \
  | --sha=* | --sh=*)
    sharedstatedir="$ac_optarg" ;;

  -site | --site | --sit)
    ac_prev=site ;;
  -site=* | --site=* | --sit=*)
    site="$ac_optarg" ;;

  -srcdir | --srcdir | --srcdi | --srcd | --src | --sr)
    ac_prev=srcdir ;;
  -srcdir=* | --srcdir=* | --srcdi=* | --srcd=* | --src=* | --sr=*)
    srcdir="$ac_optarg" ;;

  -sysconfdir | --sysconfdir | --sysconfdi | --sysconfd | --sysconf \
  | --syscon | --sysco | --sysc | --sys | --sy)
    ac_prev=sysconfdir ;;
  -sysconfdir=* | --sysconfdir=* | --sysconfdi=* | --sysconfd=* | --sysconf=* \
  | --syscon=* | --sysco=* | --sysc=* | --sys=* | --sy=*)
    sysconfdir="$ac_optarg" ;;

  -target | --target | --targe | --targ | --tar | --ta | --t)
    ac_prev=target ;;
  -target=* | --target=* | --targe=* | --targ=* | --tar=* | --ta=* | --t=*)
    target="$ac_optarg" ;;

  -v | -verbose | --verbose | --verbos | --verbo | --verb)
    verbose=yes ;;

  -version | --version | --versio | --versi | --vers)
    echo "configure generated by autoconf version AC_ACVERSION"
    exit 0 ;;


  -*) AC_MSG_ERROR([$ac_option: invalid option; use --help to show usage])
    ;;

  *)
changequote(, )dnl
    if test -n "`echo $ac_option| sed 's/[-a-z0-9.]//g'`"; then
changequote([, ])dnl
      AC_MSG_WARN($ac_option: invalid host type)
    fi
    if test "x$nonopt" != xNONE; then
      AC_MSG_ERROR(can only configure for one host and one target at a time)
    fi
    nonopt="$ac_option"
    ;;

  esac
done

if test -n "$ac_prev"; then
  AC_MSG_ERROR(missing argument to --`echo $ac_prev | sed 's/_/-/g'`)
fi
])
dnl
dnl AC_INIT_LAPS(UNIQUE-FILE-IN-SOURCE-DIR)
dnl adapted from AC_INIT autoconfig version 2.12 by Jim Edwards
dnl NOAA FSL
dnl
AC_DEFUN(AC_INIT_LAPS,
[AC_REQUIRE([AC_INIT_BINSH])dnl
AC_INIT_LAPS_NOTICE
AC_DIVERT_POP()dnl to NORMAL
AC_DIVERT_PUSH(AC_DIVERSION_INIT)dnl
AC_INIT_LAPS_PARSE_ARGS
AC_INIT_PREPARE($1)dnl
AC_DIVERT_POP()dnl to NORMAL
])

dnl AC_INIT_ARCH(ARCHNAME)
dnl
AC_DEFUN(AC_INIT_ARCH,[
#
#  Returns the arch of the machine
#
arch=$1
#
# For some reason, AIX 4.x s ARCH to blank in some cases.
if test ${ARCH}  
then
    if test -n "$ARCH"   
    then
      arch=$ARCH
    fi
fi
if test -z "$arch"  
then
    if test -x /bin/arch 
    then
       arch=`/bin/arch`
    elif test -x /bin/uname 
    then
       arch=`/bin/uname -s`
       if test "$arch" = "AIX"  
       then
             arch="rs6000"
       elif test "$arch" = "ULTRIX"  
       then
             arch=`/bin/uname -m`
             if test "$arch" = "RISC"  
             then
               arch="dec5000"
             fi
       elif test "$arch" = "HP-UX"  
       then
	     arch="hpux"
       elif test "$arch" = "IRIX64"  
       then
	     arch="IRIX"
       elif test "$arch" = "Linux"  
       then
            :
       elif test "$arch" != "IRIX"  
       then 
             arch=`/bin/uname -m`
       fi
    elif test -e /usr/local/bin/arch 
    then
         arch=`/usr/local/bin/arch`
    else
         arch="unknown"
    fi
fi


#
# Remove blanks from the value of arch (another fix for AIX 4.x)
# Note that this means that the tests below must have no blanks in them
# (this affects CRAY xxx)
#
 arch=`echo $arch | sed 's/ //g'`
 SunOSTest=`expr "$arch" : "\(....\)"`
if test "$SunOSTest" = "sun4"  
then
   arch=sun4
   Version=`/bin/uname -r`
  # In "improving" SunOS, the useful feature of "substr" was withdrawn 
  # from expr.  Can't let the users have life too easy, can we?  This 
  # means that we can't just use 
  #    MajorVersion=`expr substr $Version 1 1`
  # because it won't work on Solaris systems.  The following should work on
  # both:
   MajorVersion=`expr "$Version" : "\(.\)"`
  if test $MajorVersion -eq  5  
  then
     arch="solaris"
  fi
elif test "$arch" = "AIX" 
then
    arch="rs6000"
elif test "$arch" = "RIOS" 
then
    arch="rs6000"
elif test "$arch" = "mips" 
then
    arch="dec5000"
elif test "$arch" = "dec-5000" 
then
    arch="dec5000"
elif test "$arch" = "IP12" 
then
   arch="IRIX"
elif test "$arch" = "i386" 
then
   arch="ipsc2"
elif test "$arch" = "cray" 
then
   arch="CRAY"
elif test "$arch" = "CRAYY-MP" 
then
   arch="CRAY"
elif test "$arch" = "CRAYC90" 
then
   arch="CRAY"
elif test "$arch" = "CRAYJ90" 
then
   arch="CRAY"
elif test "$arch" = "CRAYT90" 
then
   arch="CRAY"
elif test "$arch" = "CRAYTS" 
then
   arch="CRAY"
elif test "$arch" = "CRAY-2" 
then
   arch="CRAY"
elif test "$arch" = "sun4m" 
then
   arch="sun4"
elif test "$arch" = "iris4d" 
then
   arch="IRIX"
elif test "$arch" = "next" 
then
   arch="NeXT"
elif test "$arch" = "KSR1" 
then
   arch="ksr"
elif test -x "/dev/elan" 
then
   arch="meiko"
elif test -f /usr/bin/uxpm 
then
   arch="UXPM"
fi     
#
if test -z "$arch" 
then
  AC_MSG_ERROR(Error: Couldn't guess target architecture, you must
                set an architecture type with --arch=<value>)
fi
AC_MSG_RESULT(Configuring for \"$arch\" target architecture)
CPPFLAGS="$CPPFLAGS -D$arch"

])

dnl
dnl AC_PROG_GENERIC(PROGROOT, exename, libname, incname)
dnl
AC_DEFUN(AC_PROG_GENERIC,
[

exename=$2

if test -z "[$]$1"
then
  AC_PATH_PROGS($1, $exename, ,$5)
  if test -n "[$]$1"
  then
    $1=`echo "[$]$1" | sed "s/\/bin\/$exename//"`;
  fi
fi
PROGLIB="[$]$1/lib"
PROGINC="[$]$1/include"

for ac_file in $3
do
echo "here $PROGLIB/$ac_file"
if test -f "$PROGLIB/$ac_file"
then
  AC_MSG_RESULT(Found $ac_file library in $PROGLIB)
fi
done
for ac_file in $4
do
if test -f "$PROGINC/$ac_file"
then
  AC_MSG_RESULT(Found $ac_file header in $PROGINC)
fi
done

if test -z "[$]$1"
then
  AC_MSG_ERROR(\n\nCannot find directory containing $exename in $PATH )
fi
])

---
dnl AC_PROG_GRIB2()
AC_DEFUN(AC_PROG_GRIB2,
[
dnl
dnl  Find the grib2 libraries by looking for libjasper.a, libpng.a and libz.a
dnl    in system lib locations, i.e. /usr/lib
dnl    - there must be a better way?
dnl

jtrue=0
ptrue=0
ztrue=0
if test -n "`echo $LIBS | grep ljasper`"; then jtrue=1; fi
if test -n "`echo $LIBS | grep lpng`"; then ptrue=1; fi
if test -n "`echo $LIBS | grep lz`"; then ztrue=1; fi

if test $jtrue = 1 && test $ptrue = 1 && test $ztrue = 1
then
  AC_MSG_RESULT(Found all Grib2 libraries, i.e. $LIBS)
    DEGRIBFLAGS="-DUSE_JPEG2000 -DUSE_PNG"
else
  AC_MSG_RESULT(Some Grib2 libraries -ljasper -lpng -lz were NOT FOUND ...only found >$LIBS<)
    DEGRIBFLAGS=""
fi

])

---

dnl AC_PROG_NETCDF(PATH_TO_NETCDF)
AC_DEFUN(AC_PROG_NETCDF,
[
dnl
dnl  I find the netcdf include and libraries by looking for the executables 
dnl  in the path then backing out to the corresponding lib and include dirs 
dnl  - there must be a better way?
dnl
netcdf=$1
if test -z "$netcdf" 
then
  AC_PATH_PROG(netcdf, ncdump)
  if test -n "$netcdf"
  then
    netcdf=`echo "$netcdf" | sed 's/\/bin\/ncdump//'`;
    NETCDFBIN="$netcdf/bin"
  fi
fi
NETCDFLIB="$netcdf/lib"
NETCDFINC="$netcdf/include"
if test -f $NETCDFLIB/libnetcdf.a && test -f $NETCDFINC/netcdf.h && test -f $NETCDFINC/netcdf.inc  
then
  AC_MSG_RESULT(Great Found all netcdf libraries and include files in $NETCDFLIB and $NETCDFINC)
else
  netcdf="/usr/local/netcdf"
  NETCDFBIN="$netcdf/bin"
  NETCDFLIB="$netcdf/lib"	
  NETCDFINC="$netcdf/include"
  if test -f $NETCDFLIB/libnetcdf.a && test -f $NETCDFINC/netcdf.h && test -f $NETCDFINC/netcdf.inc  
  then
    AC_MSG_RESULT(Great Found all netcdf libraries and include files in $NETCDFLIB and $NETCDFINC)
  else
    netcdf=""
  fi    
fi

if test -z "$netcdf"
then
  AC_MSG_ERROR(\n\nCannot find directory containing netcdf in $PATH 
             \nThis package is required for LAPS and must have ../lib/libnetcdf.a 
             ../include/netcdf.h and ../include/netcdf.inc
             \nIf you feel that this message is in error try setting the netcdf\n
             root using the --netcdf argument to configure)
fi


AC_SUBST(NETCDFBIN)             
AC_SUBST(NETCDFINC)
AC_SUBST(NETCDFLIB)

])

AC_DEFUN(AC_QUERY_USER_YN,
[
  while :
  do 
    read ac_query_user_yn
#    ac_query_user_yn=`grabchars -q'Answer y or n:'`
    case "$ac_query_user_yn" in
    y)  AC_MSG_RESULT(yes) 1>&2; break;;
    n)  AC_MSG_RESULT(no) 1>&2; break;;
    *)  AC_MSG_WARN($ac_query_user_yn ? Please answer y or n.) 1>&2 ;;
    esac
  done
])

AC_DEFUN(AC_QUERY_USER,
[
    read ac_query_user
    if test -z "$ac_query_user"
    then
      ac_query_user=$1
    else
      $1=$ac_query_user
    fi
])


dnl The name of shell var CACHE-ID must contain `_cv_' in order to get saved.
dnl AC_CACHE_UPDATE(CACHE-ID, COMMANDS-TO-SET-IT)
define(AC_CACHE_UPDATE,
[dnl We used to use the below line, but it fails if the 1st arg is a
dnl shell variable, so we need the eval.
dnl if test "${$1+set}" = set; then
dnl the '' avoids an AIX 4.1 sh bug ("invalid expansion").
if eval "test \"`echo '$''{'$1'+set}'`\" = set"; then
  echo $ac_n "(cached) $ac_c" 1>&AC_FD_MSG
fi
#else
  $2
#fi
])


define(LAPS_SET_PATH,
[
  AC_CACHE_UPDATE($2,[
      while :
      do
        echo "$1 ($$2)"
        AC_QUERY_USER($2)
        if test -d "$$2"
        then
          break
        else
          echo "directory $$2 not found!  Use this value anyway?"
          AC_QUERY_USER_YN
          if test "$ac_query_user_yn" = "y"
          then
            break
          fi
        fi
      done

  ])
])

dnl
dnl Configure the nest7grid.parms file for site 
dnl
AC_DEFUN(LAPS_PARMS_CONFIG,
[
  echo ""
  echo "The file $INSTALLROOT/data/static/nest7grid.parms contains runtime"  
  echo "modifiable data paths and configuration info pertaining to LAPS"
  echo "This section will take you through the file and help you to configure it."
  echo ""
  echo "Do you want to configure nest7grids.parms at this time?"
  NEST7GRID=""   
  AC_QUERY_USER_YN
  if test "$ac_query_user_yn" = "y"
  then
    echo ""
    echo "Configuring $INSTALLROOT/data/static/nest7grid.parms"
    echo ""
    NEST7GRID="data/static/nest7grid.parms"
dnl
dnl AC_CACHE_UPDATE querys the user to see if the cached value is acceptable
dnl
    AC_CACHE_UPDATE(laps_cv_stndlat,[ 
      echo "What is the standard latitude of the LAPS grid? ($laps_cv_stndlat)"
      AC_QUERY_USER(laps_cv_stndlat)
    ])
    STNDLAT=$laps_cv_stndlat    
    AC_SUBST(STNDLAT)

    AC_CACHE_UPDATE(laps_cv_centlat,[ 
      echo "What is the center latitude of the LAPS grid? ($laps_cv_centlat)"
      AC_QUERY_USER(laps_cv_centlat)
    ])
    CENTLAT=$laps_cv_centlat    
    AC_SUBST(CENTLAT)

    AC_CACHE_UPDATE(laps_cv_stndlon,[ 
      echo "What is the standard longitude of the LAPS grid? ($laps_cv_stndlon)"
      AC_QUERY_USER(laps_cv_stndlon)
    ])
    STNDLON=$laps_cv_stndlon    
    AC_SUBST(STNDLON)

    AC_CACHE_UPDATE(laps_cv_centlon,[ 
      echo "What is the center longitude of the LAPS grid? ($laps_cv_centlon)"
      AC_QUERY_USER(laps_cv_centlon)
    ])
    CENTLON=$laps_cv_centlon    
    AC_SUBST(CENTLON)


    LAPS_SET_PATH("Path to pirep data",laps_cv_path_pirep)
    PIREPPATH=$laps_cv_path_pirep
    AC_SUBST(PIREPPATH)

    LAPS_SET_PATH("Path to rass data",laps_cv_path_rass)
    RASSPATH=$laps_cv_path_rass
    AC_SUBST(RASSPATH)

    LAPS_SET_PATH("Path to profiler data",laps_cv_path_prof)
    PROFPATH=$laps_cv_path_prof
    AC_SUBST(PROFPATH)

    LAPS_SET_PATH("Path to boundary layer rass data",laps_cv_path_blprass)
    BLPRASSPATH=$laps_cv_path_blprass
    AC_SUBST(BLPRASSPATH)

    LAPS_SET_PATH("Path to boundary layer profiler data",laps_cv_path_blpprof)
    BLPPROFPATH=$laps_cv_path_blpprof
    AC_SUBST(BLPPROFPATH)

    LAPS_SET_PATH("Path to satellite cdf data",laps_cv_path_satcdf)
    SATCDF=$laps_cv_path_satcdf
    AC_SUBST(SATCDF)

    LAPS_SET_PATH("Path to RUC data",laps_cv_path_ruc)
    RUCPATH=$laps_cv_path_ruc
    AC_SUBST(RUCPATH)

    LAPS_SET_PATH("Path to NGM data",laps_cv_path_ngm)
    PIREPPATH=$laps_cv_path_ngm
    AC_SUBST(NGMPATH)


    LAPS_SET_PATH("Path to WSI 2D data",laps_cv_path_wsi2d)
    WSI2D=$laps_cv_path_wsi2d
    AC_SUBST(WSI2D)

    LAPS_SET_PATH("Path to WSI 3D data",laps_cv_path_wsi3d)
    WSI3D=$laps_cv_path_wsi3d
    AC_SUBST(WSI3D)

    LAPS_SET_PATH("Path to RAOB data",laps_cv_path_raob)
    RAOBPATH=$laps_cv_path_raob
    AC_SUBST(RAOBPATH)

    LAPS_SET_PATH("Path to QC ACARS data",laps_cv_path_qcacars)
    QCACARS=$laps_cv_path_qcacars
    AC_SUBST(QCACARS)

    LAPS_SET_PATH("Path to Satellite GVAR data",laps_cv_path_satgvar)
    SATGVAR=$laps_cv_path_satgvar
    AC_SUBST(SATGVAR)

    LAPS_SET_PATH("Path to Satellite WFO vis data",laps_cv_satwfovis)
    SATWFOVIS=$laps_cv_satwfovis
    AC_SUBST(SATWFOVIS)

    LAPS_SET_PATH("Path to Satellite WFO 39 data",laps_cv_satwfo39)
    SATWFO39=$laps_cv_satwfo39
    AC_SUBST(SATWFO39)

    LAPS_SET_PATH("Path to Satellite WFO wv data",laps_cv_satwfowv)
    SATWFOWV=$laps_cv_satwfowv
    AC_SUBST(SATWFOWV)

    LAPS_SET_PATH("Path to Satellite WFO i11 data",laps_cv_satwfoi11)
    SATWFOI11=$laps_cv_satwfoi11
    AC_SUBST(SATWFOI11)

    LAPS_SET_PATH("Path to Satellite WFO i12 data",laps_cv_satwfoi12)
    SATWFOI12=$laps_cv_satwfoi12
    AC_SUBST(SATWFOI12)

  fi
])    







