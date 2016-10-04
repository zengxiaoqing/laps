

        subroutine get_elong_a(ialt_start,ialt_end,ialt_delt  &! I
                   ,jazi_start,jazi_end,jazi_delt             &! I
                   ,minalt,maxalt,minazi,maxazi               &! I
                   ,sol_alt,sol_azi,view_alt,view_az          &! I
                   ,elong                                   )  ! O

        include 'trigd.inc'

        cosunitvectors(a1,a2,a3,b1,b2,b3) = min(max(a1*b1+a2*b2+a3*b3,-1.),+1.)
        angleunitvectors(a1,a2,a3,b1,b2,b3) = acosd(cosunitvectors(a1,a2,a3,b1,b2,b3))

        real elong(minalt:maxalt,minazi:maxazi)
        real view_alt(minalt:maxalt,minazi:maxazi)
        real view_az(minalt:maxalt,minazi:maxazi)

        do ialt = ialt_start,ialt_end,ialt_delt
  
         altray = view_alt(ialt,jazi_start)

         do jazi = jazi_start,jazi_end,jazi_delt

          altray = view_alt(ialt,jazi)
          view_altitude_deg = altray

          view_azi_deg = view_az(ialt,jazi)

          xs = cosd(sol_alt) * cosd(sol_azi)
          ys = cosd(sol_alt) * sind(sol_azi)
          zs = sind(sol_alt)

          xo = cosd(altray) * cosd(view_azi_deg)
          yo = cosd(altray) * sind(view_azi_deg)
          zo = sind(altray)

          elong(ialt,jazi) = angleunitvectors(xs,ys,zs,xo,yo,zo)    

         enddo ! jazi
        enddo ! ialt

        return
        end
