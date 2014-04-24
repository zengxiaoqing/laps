!
!  A statement function of scaling SH to 'RH' = SH/s2r(p)
!  where p is pressure in hPa.
!
!  History: Yuanfu Xie Dec. 2013.
!
	real :: s2r,ps
	s2r(ps) = 10.**(-24./log10(ps)+9.)

