
PRO dscovr_gamma

ortho_fname =  'dscovrlatest.png'

read_png,ortho_fname,earth_rgb_ortho,order=1

help,earth_rgb_ortho

earth_out = earth_rgb_ortho

idim = 2084
imax = idim-1

; Higher Gamma lowers Green (and other) output
factor = 0.95
gamma = 0.45 ; 0.7

FOR I = 0,2047 DO BEGIN
FOR J = 0,2047 DO BEGIN

;  Gamma correction
   red_old = float(earth_out[0,I,J])
   grn_old = float(earth_out[1,I,J])
   blu_old = float(earth_out[2,I,J])

   grn_gamma =  (((grn_old/256.)^gamma) * 256.) * factor          
   grn_ratio = grn_gamma / grn_old
   grn_new   = grn_gamma                

   rog_old = red_old / grn_old
   bog_old = blu_old / grn_old

   rog_new = rog_old^0.45

;  Higher power coefficient raises Blue output
;  bog_new = bog_old^0.60 ; yields 81/101/137    
;  bog_new = bog_old^0.50 ; yields 81/101/129    
   bog_new = bog_old^0.55 ; yields 81/101/129    

   red_new = grn_new * rog_new
   blu_new = grn_new * bog_new

   if(J EQ 1000)then begin
;    print,i,red_old,grn_old,blu_old,red_new,grn_new,blu_new
   endif

   earth_out[0,I,J] = (round(red_new)) < 255
   earth_out[1,I,J] = (round(grn_new)) < 255
   earth_out[2,I,J] = (round(blu_new)) < 255

;  Contrast correction
;  earth_out[IC,*,*] = ( round( (a + 50.) * 0.7  )) < 255

ENDFOR
ENDFOR

FOR I = 1,2000,50 DO BEGIN
  print,earth_rgb_ortho[1,I,1000],float(earth_rgb_ortho[1,I,1000]),earth_out[1,I,1000]
ENDFOR

outfile = 'dscovr_gamma_full.png'
write_png, outfile, earth_out[*,*,*], /ORDER   

outfile = 'dscovr_gamma_crop.png'
; increase all to enlarge, lr and tb should have the same mean value
; lcrop  = 154 ; increase to move image left
; rcrop  = 199 ; decrease to move image left
; tcrop  = 162 ; decrease to move image down
; bcrop  = 200 ; increase to move image down
imean   = 220 ; 90 
ioffset = 15 ; decrease to move image left
lcrop  = imean-ioffset 
rcrop  = imean+ioffset 

jmean  = 220 ; 90
joffset = 18 ; increase to move image up
tcrop  = jmean-joffset
bcrop  = jmean+joffset
write_png, outfile, earth_out[*,lcrop:imax-rcrop,tcrop:imax-bcrop], /ORDER   

END

