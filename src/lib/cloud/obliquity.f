
        FUNCTION Obliquity(T)

c       T is the Julian Date (TDT)

        IMPLICIT REAL*8 (A-Z)

        PI = 3.1415926535897932d0
        RPD = PI/180.d0

        Tu = (T-2451545d0)/36525d0
        O = 23.439291d0 + Tu *
     1  (-.0130042d0 + Tu * (-.00000016d0 + Tu * .000000504d0))

        Obliquity = O * RPD

        RETURN
        END

