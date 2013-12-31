! 
! Write Portable PixMap Image File (ppm, pbm, pgm) 
! with Fortran 
!
!   Charles O'Neill Oct 23, 2009
!    charles.oneill@gmail.com
!    www.caselab.okstate.edu
!
! Copyright (c) 2009 Charles O'Neill
!
! Permission is hereby granted, free of charge, to any person
! obtaining a copy of this software and associated documentation
! files (the "Software"), to deal in the Software without
! restriction, including without limitation the rights to use,
! copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the
! Software is furnished to do so, subject to the following
! conditions:
!
! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
! OTHER DEALINGS IN THE SOFTWARE.
!

module PPM
  implicit none

contains

  !--------------------------------------------------------------
  ! Portable PixMap Type 1 (Black and White)
  subroutine writeppm1Matrix(M,text)
    integer :: M(:,:)
    character(len=*) :: text
    integer :: cols,rows
    integer :: i,j
 
    ! Open File   
    open(unit=100, file=trim(text)//".pbm", status='unknown')
    
    ! Write Header and ppm file type
    write(100,'( A )') "P1"
    write(100,'( A )') "# PPM Type 1 File (generated with fortran)"
 
    ! Write Image Size
    cols = size(M,2)
    rows = size(M,1)
    write(100,'( i, 1x, i )') cols, rows
    
    ! Write Image
    do i=1,rows
      do j=1,cols
        write(100,'( i1 )', advance='no') M(i,j)
      enddo
      write(100,*) ! Endline
    enddo
  end subroutine

  !--------------------------------------------------------------
  ! Portable PixMap Type 2 (Grayscale)
  subroutine writeppm2Matrix(M,text)
    integer :: M(:,:)
    character(len=*) :: text
    integer :: cols,rows
    integer :: i,j
    integer :: maxvalue

    ! Open File   
    open(unit=100, file=trim(text)//".pgm", status='unknown')
    
    ! Write Header and ppm file type
    write(100,'( A )') "P2"
    write(100,'( A )') "# PPM Type 2 File (generated with fortran)"
 
    ! Write Image Size
    cols = size(M,2)
    rows = size(M,1)
    write(100,'( i, 1x, i )') cols, rows
    
    ! Write Maximum Value
    maxvalue = maxval(maxval(M,dim=1),dim=1)
    write(100,'( i )') maxvalue
    
    ! Write Image
    do i=1,rows
      do j=1,cols
        write(100,'( i5,1x )', advance='no') M(i,j)
      enddo
      write(100,*) ! Endline
    enddo
  end subroutine

  !--------------------------------------------------------------
  ! Portable PixMap Type 3 (RGB color)
  subroutine writeppm3Matrix(R,G,B,text)
    integer :: R(:,:),G(:,:),B(:,:)
    character(len=*) :: text
    integer :: cols,rows
    integer :: i,j
    integer :: maxvalue
 
    ! Open File   
    open(unit=100, file=trim(text)//".ppm", status='unknown')
    
    ! Write Header and ppm file type
    write(100,'( A )') "P3"
    write(100,'( A )') "# PPM Type 2 File (generated with fortran)"
 
    ! Write Image Size
    cols = size(R,2)
    rows = size(R,1)
    write(100,'( i, 1x, i )') cols, rows
    
    ! Write Maximum Value
    maxvalue = max( maxval(maxval(R,dim=1),dim=1)&
                   ,maxval(maxval(G,dim=1),dim=1)&
                   ,maxval(maxval(B,dim=1),dim=1))
    write(100,'( i )') maxvalue
    
    ! Write Image
    do i=1,rows
      do j=1,cols
        write(100,'( 3(i5,1x) )') R(i,j),G(i,j),B(i,j)
      enddo
    enddo
  end subroutine

  !--------------------------------------------------------------
  ! Test Module
  subroutine testPPM
    integer,parameter :: N = 100
    integer :: A(N,N)
    integer :: R(N,N)
    integer :: G(N,N)
    integer :: B(N,N)
    real :: AA(N,N)
    integer :: i,j

    ! Show the PixMap Format with a simple case
    open(unit=100, file="test.ppm", status='unknown')
    write(100,'( A )') "P1"
    write(100,'( A )') "# This is an example bitmap"
    write(100,'( A )') "18 10"
    write(100,'( A )') "0 0 0 0 0 0   0 0 0 0 0 0   0 0 0 0 0 0"
    write(100,'( A )') "0 0 1 1 0 0   0 0 1 1 1 0   0 1 0 0 1 0"
    write(100,'( A )') "0 1 0 0 1 0   0 1 0 0 1 0   0 1 0 0 1 0"
    write(100,'( A )') "0 1 0 0 1 0   0 1 0 0 0 0   0 1 0 0 1 0"
    write(100,'( A )') "0 1 0 0 1 0   0 0 1 0 0 0   0 1 0 0 1 0"
    write(100,'( A )') "0 1 0 0 1 0   0 0 0 1 0 0   0 1 0 0 1 0"
    write(100,'( A )') "0 1 0 0 1 0   0 0 0 0 1 0   0 1 0 0 1 0"
    write(100,'( A )') "0 1 0 0 1 0   0 1 0 0 1 0   0 1 0 0 1 0"
    write(100,'( A )') "0 0 1 1 0 0   0 1 1 1 0 0   0 0 1 1 0 0"
    write(100,'( A )') "0 0 0 0 0 0   0 0 0 0 0 0   0 0 0 0 0 0"
    close(100)
        
    ! Get a matrix of random numbers
    CALL RANDOM_SEED()
    call RANDOM_NUMBER(AA)
    
    ! Setup and write type 1
    A = AA + 0.5
    call writeppm1Matrix(A,"subtest")
    write(*,*) "writeppm1Matrix"
  
    ! Setup and write type 2
    A = AA * 256
    call writeppm2Matrix(A,"subtest2")
    write(*,*) "writeppm2Matrix"

    ! Setup and write type 3 
    do i=1,N
      do j=1,N
        R(i,j) = abs( mod(i+j,N/5) )
        G(i,j) = abs( mod(i-j,N/5) )
        B(i,j) = abs( mod(i*j,N/5) )
      enddo
    enddo
    call writeppm3Matrix(R,G,B,"subtest3")
    write(*,*) "writeppm3Matrix"
  end subroutine

end module
