
        subroutine midpoint(a1,a2,a3,b1,b2,b3,m1,m2,m3)
        implicit real*8 (A-Z)

        call normalize(a1,a2,a3,maga)
        call normalize(b1,b2,b3,magb)
        m1 = 0.5d0 * (a1+b1)
        m2 = 0.5d0 * (a2+b2)
        m3 = 0.5d0 * (a3+b3)
        call normalize(m1,m2,m3,magm)

        return
        end


        subroutine anglevectors(a1,a2,a3,b1,b2,b3,angle)
        implicit real*8 (A-Z)

        maga = mag(a1,a2,a3)
        magb = mag(b1,b2,b3)

        call crossproduct(a1,a2,a3,b1,b2,b3,c1,c2,c3)

        magc = mag(c1,c2,c3)

        sinangle = magc/(maga*magb)
        cosangle = dotproduct(a1,a2,a3,b1,b2,b3)/(maga*magb)

        angle = atan2(sinangle,cosangle)

        return
        end


        function dotproduct(a1,a2,a3,b1,b2,b3)
        implicit real*8 (A-Z)

        dotproduct  = a1 * b1 + a2 * b2 + a3 * b3

        return
        end


        function angle_vectors(a1,a2,a3,b1,b2,b3)
        implicit real*8 (A-Z)

        maga = mag(a1,a2,a3)
        magb = mag(b1,b2,b3)

        call crossproduct(a1,a2,a3,b1,b2,b3,c1,c2,c3)

        magc = mag(c1,c2,c3)

        sinangle = magc/(maga*magb)
        cosangle = dotproduct(a1,a2,a3,b1,b2,b3)/(maga*magb)

        angle_vectors = atan2(sinangle,cosangle)

        return
        end



        subroutine normalize(a1,a2,a3,maga)
        implicit real*8 (A-Z)

        maga    = sqrt(a1*a1+a2*a2+a3*a3)

        if(maga .ne. 0.)then

            magainv = 1.d0/maga

            a1 = a1 * magainv
            a2 = a2 * magainv
            a3 = a3 * magainv

        endif

        return
        end


        function mag(a1,a2,a3)
        implicit real*8 (A-Z)

        mag    = sqrt(a1*a1+a2*a2+a3*a3)

        return
        end


        subroutine xyztolatlon(x,y,z,lat,lon)
        implicit real*8 (A-Z)

        lat = atan(z/sqrt(x**2 + y**2))
        lon = atan3(y,x)

        return
        end


        subroutine latlontoxyz(lat,lon,x,y,z)
        implicit real*8 (A-Z)

        x = cos(lon)*cos(lat)
        y = sin(lon)*cos(lat)
        z = sin(lat)

        return
        end


        subroutine crossproduct(a1,a2,a3,b1,b2,b3,c1,c2,c3)
        implicit real*8 (A-Z)

        c1   = a2 * b3 - a3 * b2
        c2   = a3 * b1 - a1 * b3
        c3   = a1 * b2 - a2 * b1

        return
        end
