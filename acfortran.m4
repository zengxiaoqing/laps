dnl
dnl tests specific to fortran and it's interface to C
dnl AC_PROG_FC(additional_fortran_compilers)
dnl
AC_DEFUN(AC_PROG_FC,[
AC_CHECK_PROGS(FC,$1 f90 xlf90 pgf90 xlf f77  gf77,,$PATH)
test -z "$FC" && AC_MSG_ERROR([no acceptable fortran found in \$PATH])
cat >conftest.f <<EOF
          program main
          end
EOF
/bin/rm -f conftest.out
$FC $FFLAGS -c conftest.f > conftest.out 2>&1
if test $? != 0 ; then
    AC_MSG_RESULT(no)
    AC_MSG_WARN(Fortran compiler returned non-zero return code)
    if test -s conftest.out ; then
        cat conftest.out
    fi
elif test ! -s conftest.o ; then
    AC_MSG_RESULT(no)
    AC_MSG_WARN(Fortran compiler did not produce object file)
    if test -s conftest.out ; then
        cat conftest.out
    fi
else    
    AC_MSG_RESULT($FC seems to work)
fi
rm -f conftest* 
AC_SUBST(FC)
])

AC_DEFUN(AC_FC_MPI,[
  AC_REQUIRE([AC_PROG_FC])
  AC_MSG_CHECKING(how to compile with mpi)
  cat << EOF > mpitest.f
	program mpitest
        include 'mpif.h'
        call mpi_init
        end
EOF
/bin/rm -f mpitest.out
$FC $FFLAGS $MPIINC $MPILIB mpitest.f -o mpitest > mpitest.out 2>&1
if test $? != 0
then
  dnl
  dnl Maybe we only need to link the library
  dnl 
  /bin/rm -f mpitest.out
  $FC $FFLAGS $MPIINC $MPILIB mpitest.f -lmpi -o mpitest > mpitest.out 2>&1
  if test $? != 0
  then
    dnl Lets try to find mpirun then
    AC_PROG_GENERIC(MPI,mpirun, libmpi.a libmpi.so, mpif)
    if test -z "$MPI"
    then
      AC_MSG_ERROR(Could not identify mpi)
    else
      /bin/rm -f mpitest.out
      $FC $FFLAGS mpitest.f $MPIINC $MPILIB -I$MPI/include -L$MPI/lib -lmpi -o mpitest > mpitest.out 2>&1
      if test $? != 0
      then
        AC_MSG_ERROR(Could not compile test program with mpi in $MPI)
      else
	AC_MSG_RESULT(Adding -I$MPI/include to MPIINC and  -L$MPI/lib -lmpi to MPILIB)
        MPILIB="$MPILIB  -L$MPI/lib -lmpi"
        MPIINC="$MPIINC -I$MPI/include"
      fi
    fi      
  else
    MPILIB="$MPILIB -lmpi"
    AC_MSG_RESULT(Adding -lmpi to MPILIB)
  fi
else
  AC_MSG_RESULT(No special requirments)
fi
AC_SUBST(MPIINC)
AC_SUBST(MPILIB)
])

AC_DEFUN(AC_FC_FUNC_ALLOC,[
  AC_REQUIRE([AC_PROG_FC])dnl
  AC_MSG_CHECKING(For allocate feature in fortran)
  cat << EOF > testalloc.f
       program testalloc
       real,allocatable:: a(:,:)
       integer i
       
       read(5,*) i
       allocate(A(i,i))
       deallocate(a)
       print*, 'allocate works'
       end
EOF
/bin/rm -f testalloc.out
$FC $FFLAGS testalloc.f -o testalloc > testalloc.out 2>&1
if test $? != 0
then
  AC_MSG_RESULT(no)
  HAVE_FC_ALLOC=0
else
  AC_MSG_RESULT(yes)
  HAVE_FC_ALLOC=1
fi
rm -f testalloc.*
])


AC_DEFUN(AC_FC_BYTE_UNSIGNED,[
  AC_REQUIRE([AC_PROG_FC])dnl
  AC_MSG_CHECKING(to see if $FC recognizes type byte as unsigned)
  rm -f testbyte*
  cat << EOF > testbyte.f
       program testbyte
       byte a
       data a/255/
       print *, a
       end
EOF
$FC $FFLAGS testbyte.f -o testbyte > testbyte.out 2>&1
if test $? != 0
then
  AC_MSG_RESULT(no)
  HAVE_FC_BYTE=0
else
  AC_MSG_RESULT(yes)
  FC_BYTE_MISSING=255
  HAVE_FC_BYTE=1
fi
rm -f testbyte*
])


AC_DEFUN(AC_FC_BYTE_SIGNED,[
  AC_REQUIRE([AC_PROG_FC])dnl
  AC_MSG_CHECKING(to see if $FC recognizes type byte as signed)
  /bin/rm -f testbyte*
  cat << EOF > testbyte.f
       program testbyte
       byte a
       data a/-1/
       print *, a
       end
EOF
  /bin/rm -f testbyte.out
  $FC $FFLAGS testbyte.f -o testbyte > testbyte.out 2>&1
  if test $? != 0	
  then
    AC_MSG_RESULT(no)
    HAVE_FC_BYTE=0
  else
    AC_MSG_RESULT(yes)
    FC_BYTE_MISSING=-1
    HAVE_FC_BYTE=1
  fi

])



AC_DEFUN(AC_FC_INT1,[
  AC_REQUIRE([AC_PROG_FC])dnl
  AC_MSG_CHECKING(to see if $FC recognizes type integer*1)

  /bin/rm -f testbyte*
  cat << EOF > testbyte.f
       program testbyte
       integer*1 a
       data a/-1/
       print *, a
       end
EOF
  /bin/rm -f testbyte.out
  $FC $FFLAGS testbyte.f -o testbyte > testbyte.out 2>&1
  if test $? != 0	
  then
    AC_MSG_RESULT(no)
    HAVE_FC_INT1=0
  else
    AC_MSG_RESULT(yes)
    FC_BYTE_MISSING=-1
    HAVE_FC_INT1=1
  fi

])


AC_DEFUN(AC_F90_FUNC_TRIGD,[
  AC_REQUIRE([AC_PROG_FC])dnl
  AC_MSG_CHECKING(For intrinsic degree based trig functions in $FC)
	cat <<EOF > testtrigd.f
	  program t_trigd
	  real a
	  a = sind(45.)
	  print *, a
	end
EOF
  /bin/rm -f testtrigd.out
  $FC $FFLAGS testtrigd.f -o testtrigd > testtrigd.out 2>&1
  if test $? != 0
  then
    AC_MSG_RESULT(no)
    USE_TRIGD=1
    FC_USE_TRIGD='      use trigd'
  else
    AC_MSG_RESULT(yes)
    USE_TRIGD=0
    FC_USE_TRIGD='C      use trigd'
  fi
  rm -f testtrigd*
])

AC_DEFUN(AC_FC_FUNC_NCARGFC,[
  AC_REQUIRE([AC_PROG_FC])dnl
  AC_MSG_CHECKING(how to compile with ncar graphics)
  AC_PATH_PROGS(NCARGFC,ncargf77)
  if test "$NCARGFC"  
  then
        cat <<EOF > testncarg_sub.f
          subroutine tncarg_sub
          real xplot, py1
          CALL FRSTPT(XPLOT,PY1)
          return
          end
EOF
	cat <<EOF > testncarg.f
	  program t_ncarg
	  real a
	  call opngks
          call tncarg_sub
	  call clsgks
	end
EOF
  /bin/rm -f testncarg.out
  $FC $FFLAGS -c testncarg_sub.f 1>  testncarg_sub.out 2>&1
  if test $? = 0 
  then
    $NCARGFC $FFLAGS testncarg.f testncarg_sub.o -o testncarg > testncarg.out 2>&1
  fi
  if test $? = 0
  then
    AC_MSG_RESULT("It appears that ncarg will work as is.")
  else
    AC_MSG_RESULT("\n $NCARGFC will not compile a test prog as is.  Try running ncargf90.pl before you compile")
    NCARGFC="ncargf90"
  fi
  rm -f testncarg*

  fi
  if test -z "$NCARGFC"
  then
    AC_MSG_WARN(Could not find NCAR Graphics. If you wish to use this package you must have the environment variable NCARG_ROOT set and the ncargf77 program in your path)
    NCARGFC="$FC"
  fi
])


AC_DEFUN(AC_FC_FUNC_GETENV,[
  AC_REQUIRE([AC_PROG_FC])dnl
  AC_MSG_CHECKING(For a getenv function in fortran)
  cat << EOF > testgetenv.f
       program t_getenv
       character*8 name
       character*80 response
       data name/'HOME'/
       call getenv(name,response)
       print *,response
       end
EOF
/bin/rm -f testgenv.out
$FC $FFLAGS testgetenv.f -o genv > testgenv.out 2>&1
if test $? != 0
then
  AC_MSG_RESULT(no)
  HAVE_FC_GETENV=0
else
  AC_MSG_RESULT(yes)
  HAVE_FC_GETENV=1
  THOME=`genv`
  if test ! "$THOME"="$HOME"
  then
    AC_MSG_ERROR(Fortran getenv does not seem to return the proper value for 
                 the environment variable HOME env=$HOME $FC=$THOME)
  fi
fi
rm -f testgetenv.f genv testgenv.out  
])

dnl Fortran runtime for Fortran/C linking
dnl On suns, try
dnl FC_LIB          =/usr/local/lang/SC2.0.1/libM77.a \ 
dnl              /usr/local/lang/SC2.0.1/libF77.a -lm \
dnl              /usr/local/lang/SC2.0.1/libm.a \
dnl              /usr/local/lang/SC2.0.1/libansi.a
dnl
dnl AIX requires -bI:/usr/lpp/xlf/lib/lowsys.exp
dnl ------------------------------------------------------------------------
dnl
dnl Get the format of Fortran names.  Uses F77, FFLAGS, and sets WDEF.
dnl If the test fails, sets NOF77 to 1, HAS_FORTRAN to 0
dnl
define(AC_FC_NAMES,[
  AC_REQUIRE([AC_PROG_FC])dnl
   # Check for strange behavior of Fortran.  For example, some FreeBSD
   # systems use f2c to implement f77, and the version of f2c that they 
   # use generates TWO (!!!) trailing underscores
   #
   # Eventually, we want to be able to override the choices here and
   # force a particular form.  This is particularly useful in systems
   # where a Fortran compiler option is used to force a particular
   # external name format (rs6000 xlf, for example).
   cat > confftest.f <<EOF
       subroutine iNiT_fop( a )
       integer a
       a = 1
       return
       end
EOF
   $FC $FFLAGS -c confftest.f > /dev/null 2>&1
   if test ! -s confftest.o ; then
        AC_MSG_ERROR(Unable to test Fortran compiler:
        Compiling a test program failed to produce an object file.)
   elif test -n "$FORTRANNAMES" ; then
	WDEF="-D$FORTRANNAMES"
   else
    # We have to be careful here, since the name may occur in several
    # forms.  We try to handle this by testing for several forms
    # directly.
    if test $arch_CRAY ; then
     # Cray doesn't accept -a ...
     nameform1=`strings confftest.o | grep init_fop_  | head -1`
     nameform2=`strings confftest.o | grep INIT_FOP   | head -1`
     nameform3=`strings confftest.o | grep init_fop   | head -1`
     nameform4=`strings confftest.o | grep init_fop__ | head -1`
    else
     nameform1=`strings -a confftest.o | grep init_fop_  | head -1`
     nameform2=`strings -a confftest.o | grep INIT_FOP   | head -1`
     nameform3=`strings -a confftest.o | grep init_fop   | head -1`
     nameform4=`strings -a confftest.o | grep init_fop__ | head -1`
    fi
    /bin/rm -f confftest.f confftest.o
    if test -n "$nameform4" ; then
	AC_MSG_RESULT(Fortran externals are lower case and have 1 or 2 trailing underscores)
	WDEF=-DFORTRANDOUBLEUNDERSCORE
    elif test -n "$nameform1" ; then
        AC_MSG_RESULT(Fortran externals have a trailing underscore and are lowercase)
	WDEF=-DFORTRANUNDERSCORE
    elif test -n "$nameform2" ; then
	AC_MSG_RESULT(Fortran externals are uppercase)
	WDEF=-DFORTRANCAPS 
    elif test -n "$nameform3" ; then
	AC_MSG_RESULT(Fortran externals are lower case)
	WDEF=-DFORTRANNOUNDERSCORE 
    else
	AC_MSG_ERROR(Unable to determine the form of Fortran external names
          Make sure that the compiler $FC can be run on this system)
    fi
    fi])dnl


AC_DEFUN(AC_FC_FUNC_IMPLICIT_ALLOC,[
  AC_REQUIRE([AC_PROG_FC])dnl
# 
# Check for implicit dynamic memory allocation support in fortran
#
# This is unnessaccerily complicated in an attempt to fool optimizers 
#
   cat > memtest.f <<EOF
       program memtest
       integer n1, n2
       n1=200
       n2=300
       call mem_test(n1,n2)
       
       end
       subroutine mem_test( n1, n2 )
       integer n1,n2, i, j
       real f(n1,n2)
       do j=1,n2
         do i=1,n1
           f(i,j) = i+j*n1
         enddo
       enddo
       print *,f
       return
       end
EOF
$FC $FFLAGS  memtest.f -o memtest > /dev/null 2>&1
if test ! -x memtest ; then
  AC_MSG_RESULT(Test for implicit memory allocation in fortran did not compile
                   assuming implicit allocation is not supported)
  DYNAMIC=0
else
  memtest > /dev/null 2>&1
  if test $? != 0 
  then
    AC_MSG_RESULT(Test for implicit memory allocation in fortran failed
                   assuming implicit allocation is not supported)
    DYNAMIC=0
  else
    AC_MSG_RESULT(Good: $FC supports implicit memory allocation)
    DYNAMIC=1
  fi
fi
rm -f memtest.* memtest
])
dnl
dnl How to use the preprocessor with fortran
dnl 
AC_DEFUN(AC_FC_CPP,[
  AC_REQUIRE([AC_PROG_FC])dnl
  AC_REQUIRE([AC_PROG_CC])dnl
  rm -f cpptest.* cpptest
  cat > cpptest.F << EOF
       program cpptest
       integer a
       a=-1
#if defined(CPPTEST)
       a=0
#endif
       print *,a
       end
EOF
dnl
dnl We want to make sure we get a clean fortran output from cpp
dnl

CPP="$FC $FFLAGS"
$CPP $CPPFLAGS -DCPPTEST cpptest.F -o cpptest 1>/dev/null 2>&1
if test -x cpptest
then
  fc_tmp=`cpptest`
  USECPP=''
else
  AC_MSG_WARN($FC $FFLAGS does not seem to handle cpp flags)
  CPP="$CC -P"
  $CPP $CPPFLAGS -DCPPTEST cpptest.F 1>/dev/null 2>&1
  if test -s cpptest.i
  then
    mv cpptest.i cpptest.f
    $FC $FFLAGS cpptest.f -o cpptest 1>/dev/null 2>&1
    if test -x cpptest
    then
      fc_tmp=`cpptest`
    fi
    AC_MSG_RESULT($CPP $CPPFLAGS seems to work)  
    USECPP='USECPP=1'
  fi

fi


AC_SUBST(USECPP)
rm -f cpptest.* cpptest
])

AC_DEFUN(AC_FC_INC,[
dnl test -I command line option
  AC_REQUIRE([AC_PROG_FC])dnl
  cat > /tmp/inctest.inc << EOF
       integer a
       data a/1/
EOF
  cat > inctest.f << EOF
       program inctest
       include 'inctest.inc'
       print *,a
       end
EOF
dnl
dnl  
dnl
$FC $FFLAGS -I/tmp inctest.f -o inctest 1>/dev/null 2>&1
if test ! -x inctest ; then
  AC_MSG_RESULT(It appears that -I cannot be used with Fortran include
	        statement. Will add links for include files)

  ac_fc_inc_val=0
else
  inctest > /dev/null 2>&1
  if test $? != 0 
  then
    AC_MSG_RESULT(It appears that -I cannot be used with Fortran include
	        statement. Will add links for include files)
    ac_fc_inc_val=0
  else
    AC_MSG_RESULT(Good: $FC supports -I option for fortran include statement)
    ac_fc_inc_val=1
  fi
fi
rm -f inctest.f inctest /tmp/inctest.inc

])










