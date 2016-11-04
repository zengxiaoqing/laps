pro compare_allsky

print,' Start IDL Procedure compare_allsky (for polar)...'

rpd = 3.14159265 / 180.

filetime = '2013-07-24T19h22m47.858s-raw.jpg'

filetime = GETENV('LAPS_ASCIITIME')
site = GETENV('SITE')
file_model_polar = GETENV('FILE_MODEL_POLAR')
file_camera_polar = GETENV('FILE_CAMERA_POLAR')
a9time = GETENV('LAPS_A9TIME')
solar_alt = GETENV('SOLAR_ALT')
solar_az  = GETENV('SOLAR_AZ')
lapsdataroot = GETENV('LAPS_DATA_ROOT')

vmiss = -99.0

print, 'site from env is: ',site
print, 'solar alt/az is: ',solar_alt,'  ',solar_az

;rilaps_cen = 2.0
;rjlaps_cen = 3.0
;azi_r = atan(-rilaps_cen,+rjlaps_cen)
;print,'test atan ',rilaps_cen,rjlaps_cen,azi_r

print,' Read camera polar image from ',file_camera_polar
read_png,file_camera_polar,camera_polar
help,camera_polar

print,' Read model polar image from ',file_model_polar
read_png,file_model_polar,model_polar
help,model_polar

print,' lapsdataroot = ',lapsdataroot

idim = 511L     
jdim = 511L 

imax = idim-1
jmax = jdim-1

; Dimension mask arrays
cont_table_mask = bytarr(3,idim,jdim)
cont_table_mask[*,*,*] = 0 ; initialize black (clear sky)

camera_mask = bytarr(3,idim,jdim)
camera_mask[0,*,*] = 95  ; initialize blue (clear sky)
camera_mask[1,*,*] = 95  ; initialize blue (clear sky)
camera_mask[2,*,*] = 255 ; initialize blue (clear sky)

; Define alt/azi arrays for polar image

alt_d_a = fltarr(idim,jdim)
azi_d_a = fltarr(idim,jdim)

icont      = lonarr(2,2)
icont_rand = lonarr(2,2)

ncloud_camera = 0L
ncloud_model = 0L
npotl = 0L

icont[*,*] = 0L
icont_rand[*,*] = 0L

terralt = fltarr(13)

if(site EQ 'las')then begin
;               0          90           180            270            360
    terralt = [10., 5., 5., 5., 5.,  7.,  7.,  5.,  5.,  5.,  5.,  5., 10.]
    blue_camera_thresh = 20. ; 1.40
    thresh_pix = 0

endif else if(site EQ 'dsrc')then begin
;               0          90           180            270            360
    terralt = [ 5., 5., 5., 5., 5.,  7.,  7., 12., 13., 12., 10.,  6.,  5.]
    if(solar_alt GT 20)then begin
        blue_camera_thresh = 1.30
    endif else begin
        blue_camera_thresh = 1.15
    endelse
    thresh_pix = 0

endif else if(site EQ 'mtevans')then begin
;               0          90           180            270            360
    terralt = [ 2., 2., 2., 2., 2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.]
    blue_camera_thresh = 20. ; 1.30
    if(solar_alt GT -4.0)then begin ; day  
        thresh_pix = 0
    endif else begin
        thresh_pix = 999
    endelse
endif    
print, 'blue_camera_thresh is:',blue_camera_thresh
print, 'thresh_pix is:        ',thresh_pix           

;azi_r = atan(-rilaps_cen,+rjlaps_cen)
;print,'test atan ',rilaps_cen,rjlaps_cen,azi_r

for i = 0,imax    do begin
for j = jmax,0,-1 do begin
    ip = (i+1) < imax
    im = (i-1) > 0
    jp = (j+1) < jmax
    jm = (j-1) > 0

    rilaps_cen = float(i - 256.)
    rjlaps_cen = float(j - 256.)
    azi_r = atan(-rilaps_cen,+rjlaps_cen)

;    if (I EQ 200) then begin
;        azi_r = atan(-rilaps_cen,+rjlaps_cen)
;        print,'test atan ',rilaps_cen,rjlaps_cen,azi_r
;    endif

    azi_d = azi_r / rpd
    if(azi_d LT 0.)then begin
        azi_d = azi_d + 360.
    endif

    sky_radii_laps = sqrt(rilaps_cen^2 + rjlaps_cen^2) / 256.

    alt_d = 90. - (sky_radii_laps * 90.) 

    elg_d = acos( sin(alt_d*rpd) * sin(solar_alt*rpd) + cos(alt_d*rpd) * cos(solar_alt*rpd) * cos((azi_d-solar_az)*rpd) ) / rpd

    itn = round(azi_d / 30.)

    terrain_alt = terralt[itn]

;   Center point
    max_pix  =  (camera_polar[0,i,j] > camera_polar[1,i,j]) > camera_polar[2,i,j]                    

;   Add Surrounding points
    max_pix2 = ((camera_polar[0,ip,j] > camera_polar[1,ip,j]) > camera_polar[2,ip,j]) > max_pix
    max_pix2 = ((camera_polar[0,i,jp] > camera_polar[1,i,jp]) > camera_polar[2,i,jp]) > max_pix2

;   Almucantar brightness
    if(J EQ 256 AND I EQ (I/10)*10)then begin ; debug output
;   if(abs(alt_d - solar_alt) LT 0.2)then begin
        print,'maxpix2 ',alt_d,azi_d,elg_d,max_pix2
    endif

    if(elg_d GE 15 AND alt_d GE terrain_alt AND alt_d GE 15.0 AND max_pix2 LT 200 AND max_pix GT thresh_pix)then begin

;       color_convert,camera_polar[0,i,j],camera_polar[1,i,j],camera_polar[2,i,j],huec,satc,valuec,/RGB_HSV
;       color_convert, model_polar[0,i,j], model_polar[1,i,j], model_polar[2,i,j],huem,satm,valuem,/RGB_HSV

        numerator_camera = float(camera_polar[2,i,j])
        numerator_model = float(model_polar[2,i,j])

;       denom_camera = ( ( float(camera_polar[0,i,j] + camera_polar[1,i,j]) ) * 0.5 )
;       denom_camera = (float(camera_polar[0,i,j]) + float(camera_polar[1,i,j])) * 0.5
        denom_camera =                               float(camera_polar[1,i,j])     

;       denom_model = ( ( float(model_polar[0,i,j]  + model_polar[1,i,j] ) ) * 0.5 )
        denom_model  = (float(model_polar[0,i,j]) + float(model_polar[1,i,j])) * 0.5

        blue_camera = numerator_camera / denom_camera
;       blue_camera = float(camera_polar[2,i,j]) - float(camera_polar[1,i,j])
        blue_model  = numerator_model  / denom_model

        if(solar_alt GT -4.0)then begin ; day  
            if(blue_camera LT blue_camera_thresh)then begin
                icloud_camera = 1
            endif else begin
                icloud_camera = 0
            endelse

            if(blue_model LT 1.30)then begin
                icloud_model = 1
            endif else begin
                icloud_model = 0
            endelse

        endif else begin               ; night
            max_model = (model_polar[0,i,j] > model_polar[1,i,j]) > model_polar[2,i,j]
            max_camera = max_pix

            if(max_camera GT 128)then begin
                icloud_camera = 1
            endif else begin
                icloud_camera = 0
            endelse

            if(max_model GT 128)then begin
                icloud_model = 1
            endif else begin
                icloud_model = 0
            endelse

        endelse
        
;       Note that correct clear is black for the contingency table (as initialized)

        if(icloud_model GT icloud_camera)then begin                       ; overforecast is red
            cont_table_mask[0,i,j] = 255 
            cont_table_mask[1,i,j] =  95 
            cont_table_mask[2,i,j] =  95 

            camera_mask[0,i,j] =  95                                      ; camera has clear 
            camera_mask[1,i,j] =  95 
            camera_mask[2,i,j] = 255 
        endif

        if(icloud_model LT icloud_camera)then begin                       ; underforecast is blue
            cont_table_mask[0,i,j] =  95 
            cont_table_mask[1,i,j] =  95 
            cont_table_mask[2,i,j] = 255 

            camera_mask[0,i,j] = 220                                      ; camera has cloud
            camera_mask[1,i,j] = 220 
            camera_mask[2,i,j] = 220 
        endif

        if(icloud_camera * icloud_model GT 0)then begin
            cont_table_mask[*,i,j] = 220                                  ; correct hit is white

            camera_mask[0,i,j] = 220                                      ; camera has cloud
            camera_mask[1,i,j] = 220 
            camera_mask[2,i,j] = 220 
        endif

;       Observed is cloud   (Blue): 220,255
;       Observed is clear   (Blue): 0, 95
;       Observed is unknown (Blue): 60

        if(J EQ 256 AND I EQ (I/10)*10)then begin ; debug output
;       if(abs(alt_d - solar_alt) LT 0.2)then begin
;           print,'model polar RG',float(model_polar[0,i,j]),float(model_polar[1,i,j])
;           print,'camera polar RG',float(camera_polar[0,i,j]),float(camera_polar[1,i,j])
;           print,'num/denom',numerator_camera,numerator_model,denom_camera,denom_model
            print,FORMAT = '("i,j,alt,azi,elg,teralt,camera/model RGB,camera/model blue icloud",2I4,"  ",4F6.1,"  ",3I4,"  ",3I4,"  ",2F8.3,"  ",2I3)' $
                             ,i,j,alt_d,azi_d,elg_d,terrain_alt,camera_polar[0:2,i,j],model_polar[0:2,i,j],blue_camera,blue_model,icloud_camera,icloud_model
        endif

;       Note that when clouds are present the contingency table has a zero value
        icont[1-icloud_camera,1-icloud_model] = icont[1-icloud_camera,1-icloud_model] + 1

        npotl = npotl + 1

        ncloud_camera = ncloud_camera + icloud_camera
        ncloud_model  = ncloud_model  + icloud_model

    endif else begin ; unknown point
        imask_red = 60
        imask_grn = 80
        imask_blu = 60
        cont_table_mask[0,i,j] = imask_red                                      ; masked out area is greenish gray
        cont_table_mask[1,i,j] = imask_grn                                      ; masked out area is greenish gray
        cont_table_mask[2,i,j] = imask_blu                                      ; masked out area is greenish gray

        camera_mask[0,i,j] = imask_red                                          ; masked out area is greenish gray
        camera_mask[1,i,j] = imask_grn                                          ; masked out area is greenish gray
        camera_mask[2,i,j] = imask_blu                                          ; masked out area is greenish gray

    endelse


endfor
endfor

print,' ncloud_camera / ncloud_model = ',ncloud_camera,ncloud_model
print,' npotl = ',npotl

frac_cloud_camera = float(ncloud_camera) / float(npotl)
frac_cloud_model = float(ncloud_model) / float(npotl)

print,' frac_cloud_camera / frac_cloud_model = ',frac_cloud_camera,frac_cloud_model

print,' icont = ',icont

print,'                        Obs'
print,'                    Y          N'
print,'  Fcst Y  ',icont[0,0],icont[1,0]
print,'  Fcst N  ',icont[0,1],icont[1,1]

hits              = icont[0,0]
misses            = icont[0,1]
false_alarms      = icont[1,0]
correct_negatives = icont[1,1]

total = hits + misses + false_alarms + correct_negatives

; Accuracy
if(total GT 0)then begin
    frac_negatives = float(correct_negatives) / float(total)
    frac_coverage = 1.0 - frac_negatives

    frac_obs  = float(hits + misses)       / float(total)
    frac_fcst = float(hits + false_alarms) / float(total)

    accuracy = float(hits + correct_negatives) / float(total)

    frac_random = frac_fcst

;                           Cloudy Part                        Clear Part
;                      LAPS cloudy  Actual Cloudy    LAPS Clear             Actual Clear
    accuracy_random = (frac_random * frac_obs) + ( (1.0 - frac_random) * (1.0 - frac_obs) )
endif else begin
    frac_obs = vmiss
    frac_fcst = vmiss
    accuracy = vmiss
    accuracy_random = vmiss
endelse

; ETS
if(total GT 0)then begin
    hits_random = ((float(hits)+float(misses)) * (float(hits)+float(false_alarms))) / float(total)

    ets_denom = float(hits) + float(misses) + float(false_alarms) - hits_random

    if(ets_denom GT 0.)then begin
        ets = (float(hits) - hits_random) / ets_denom
    endif else begin
        ets = -999.9
    endelse
endif else begin
    ets = vmiss  
endelse

print,' Percent observed cloudy is ',frac_obs*100.
print,' Percent analyzed cloudy is ',frac_fcst*100.
print,' Percent correct is         ',accuracy*100.

print,' Percent correct random is  ',accuracy_random*100.

;   Write out contingency table mask image
dirout = lapsdataroot + '/' + 'lapsprd/verif/allsky/stats'
maskfile = dirout + '/' + 'verif_allsky_mask' + '.' + site + '.' + a9time + '.png'
print,' maskfile is  ',maskfile
WRITE_PNG,maskfile,cont_table_mask

;   Write out camera mask image
dirout = lapsdataroot + '/' + 'lapsprd/verif/allsky/stats'
maskfile = dirout + '/' + 'camera_allsky_mask' + '.' + site + '.' + a9time + '.png'
print,' maskfile is  ',maskfile
WRITE_PNG,maskfile,camera_mask

;   Compare to random analysis having either 50% or the actual percentage of clouds
; dirout = '/data/fab/dlaps/projects/roc/hires2/log'
statsfile = 'verif_allsky_anal' + '.' + site + '.' + a9time
print,' full statsfile is  ',dirout + '/' + statsfile

openw, 1, dirout + '/' + statsfile

printf,1,'       hits       misses     false_alarms correct_negatives total'
printf,1,        hits,      misses,    false_alarms,correct_negatives,total

printf,1,'     frac_obs     frac_fcst     accuracy   accuracy_random'
printf,1,      frac_obs,    frac_fcst,    accuracy,  accuracy_random

bias = (float(hits) + float(false_alarms)) / (float(hits) + float(misses))

printf,1,FORMAT = '("Bias: ",F7.3,"    ETS: ",F7.3)',bias,ets
printf,1,FORMAT = '(A16,4f10.4," gnuplot")',filetime,frac_obs,frac_fcst,accuracy,accuracy_random
printf,1,FORMAT = '("    solar_alt: ",F7.3)',solar_alt
printf,1,FORMAT = '("    solar_az: ",F7.3)',solar_az

ijdim = idim*jdim
print,' dimension of correlation array is  ',ijdim
for ic = 0,2 do begin
    slow = 0L & shigh = 0L & scr = -1L                  
    model_rebin = intarr(ijdim)
    help,model_rebin

    camera_rebin = intarr(ijdim)
    help,camera_rebin

    for ii = 0,imax do begin
      for jj = 0,jmax do begin
        if(cont_table_mask[1,ii,jj] NE imask_grn)then begin ; not a masked point
          scr = scr + 1
          model_rebin[scr]   = model_polar[ic,ii,jj]             
          camera_rebin[scr]  = camera_polar[ic,ii,jj]             
          if(jj EQ 256)then begin
            print,ic,ii,scr,model_rebin(scr),camera_rebin(scr)
          endif
        endif
      endfor
    endfor
    print,    'total scr is ',scr
    print,    'correlation of color ',ic,CORRELATE(model_rebin[0:scr],camera_rebin[0:scr])
    printf,1, 'correlation of color ',ic,CORRELATE(model_rebin[0:scr],camera_rebin[0:scr])
endfor

close,1

print,' compare_allsky.pro completed...'

exit

end
