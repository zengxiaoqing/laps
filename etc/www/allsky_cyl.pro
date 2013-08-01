pro allsky_cyl

; Used by the LAPS on-the-fly page

print,' Start IDL Procedure allsky_cyl...'

npts = 10000; idim * jdim

MODE_POLAR = 0
MODE_CYL = 1

IF (MODE_POLAR EQ 1) THEN BEGIN
  idim = 511     
  jdim = 511 
  allsky_polar_rgb = bytarr(3,idim,jdim)
  allsky_polar_rgb_out = bytarr(3,idim,jdim)
ENDIF ELSE BEGIN ; MODE_CYL EQ 1
  idim = 361     
; Obtain jdim from ENV
  jdim = GETENV('ALLSKY_JDIM')
  print,'Preliminary jdim from ALLSKY_JDIM environment variable is ',jdim
  jdim = 91; 92  
  allsky_polar_rgb = bytarr(3,jdim,idim)
  allsky_polar_rgb_out = bytarr(3,idim,jdim)
ENDELSE

imax = idim-1
jmax = jdim-1

for LOOP = 17,18 do begin
  
  IF (MODE_POLAR EQ 1) THEN BEGIN
    if(LOOP EQ 17)then begin & outname='allsky_polar_001.png' & allsky_rgb_name = 'allsky_rgb_polar.001' & endif
    if(LOOP EQ 18)then begin & outname='allsky_polar_002.png' & endif
  ENDIF ELSE BEGIN ; MODE_CYL EQ 1
    if(LOOP EQ 17)then begin & outname='allsky_polar_001.png' & allsky_rgb_name = 'allsky_rgb_cyl.001' & endif
    if(LOOP EQ 18)then begin & outname='allsky_polar_002.png' & endif
  ENDELSE

  if(LOOP LE 16)then begin
    istep = 30 & jstep = 10
  endif else begin
    istep = 50 & jstep = 30
  endelse

  IF (MODE_POLAR EQ 1 OR MODE_CYL EQ 1) THEN BEGIN
    print,' Reading ',allsky_rgb_name
    OPENR,3,allsky_rgb_name
    READF,3,allsky_polar_rgb
    CLOSE,3 
    help,allsky_polar_rgb

    for J = 0,jmax do begin
      for I = 0,imax do begin
        allsky_polar_rgb_out[*,I,J] = allsky_polar_rgb[*,J,I]
      endfor
    endfor
    print,' Writing ',outname
    write_png,outname,allsky_polar_rgb_out                
  ENDIF ELSE BEGIN ; MODE_CYL EQ 1
  ENDELSE

endfor

end
