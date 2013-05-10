

subroutine interp_1d(xvals,yvals,nvals,x,y,r_missing_data,idebug)

real xvals(nvals)
real yvals(nvals)
real x,y

if(idebug .eq. 1)then
    write(13,*)' x,nvals,xvals ',x,nvals,xvals
endif

y = r_missing_data

do i = 1,nvals-1
    if(xvals(i) .ge. x .and. xvals(i+1) .le. x)then
        frac = (x-xvals(i)) / (xvals(i+1)-xvals(i))
        y = yvals(i) * (1.-frac) + yvals(i+1)*frac
    endif
enddo

return
end
