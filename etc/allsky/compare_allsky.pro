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

print,' Read camera polar image from ',file_camera_polar
read_png,file_camera_polar,camera_polar
help,camera_polar

print,' Read model polar image from ',file_model_polar
read_png,file_model_polar,model_polar
help,model_polar

print,' lapsdataroot = ',lapsdataroot

idim = 511     
jdim = 511 

imax = idim-1
jmax = jdim-1

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
    terralt = [ 5., 5., 5., 5., 5.,  7.,  7., 10., 13., 10., 10.,  6.,  5.]
    if(solar_alt GT 20)then begin
        blue_camera_thresh = 20. ; 1.30
    endif else begin
        blue_camera_thresh = 20. ; 1.15
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

for i = 0,imax    do begin
for j = jmax,0,-1 do begin
    rilaps_cen = i - 256.
    rjlaps_cen = j - 256.
    azi_r = atan(rilaps_cen,rjlaps_cen)
    azi_d = azi_r / rpd
    if(azi_d lt 0.)then begin
        azi_d = azi_d + 360.
    endif

    sky_radii_laps = sqrt(rilaps_cen^2 + rjlaps_cen^2) / 256.

    alt_d = 90. - (sky_radii_laps * 90.) 

    itn = round(azi_d / 30.)

    terrain_alt = terralt[itn]

    max_pix = (camera_polar[0,i,j] > camera_polar[1,i,j]) > camera_polar[2,i,j]                    

    if(alt_d GE terrain_alt AND max_pix LT 253 AND max_pix GT thresh_pix)then begin

;       color_convert,camera_polar[0,i,j],camera_polar[1,i,j],camera_polar[2,i,j],huec,satc,valuec,/RGB_HSV
;       color_convert, model_polar[0,i,j], model_polar[1,i,j], model_polar[2,i,j],huem,satm,valuem,/RGB_HSV

        numerator_camera = float(camera_polar[2,i,j])
        numerator_model = float(model_polar[2,i,j])

;       denom_camera = ( ( float(camera_polar[0,i,j] + camera_polar[1,i,j]) ) * 0.5 )
;       denom_camera = (float(camera_polar[0,i,j]) + float(camera_polar[1,i,j])) * 0.5
        denom_camera =                               float(camera_polar[1,i,j])     

;       denom_model = ( ( float(model_polar[0,i,j]  + model_polar[1,i,j] ) ) * 0.5 )
        denom_model  = (float(model_polar[0,i,j]) + float(model_polar[1,i,j])) * 0.5

;       blue_camera = numerator_camera / denom_camera
        blue_camera = float(camera_polar[2,i,j]) - float(camera_polar[1,i,j])
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

        if(i EQ 256)then begin
;           print,'model polar RG',float(model_polar[0,i,j]),float(model_polar[1,i,j])
;           print,'camera polar RG',float(camera_polar[0,i,j]),float(camera_polar[1,i,j])
;           print,'num/denom',numerator_camera,numerator_model,denom_camera,denom_model
            print,FORMAT = '("i,j,alt,azi,teralt,camera/model RGB,camera/model blue icloud",2I4,"  ",3F6.1,"  ",3I4,"  ",3I4,"  ",2F8.3,"  ",2I3)' $
                             ,i,j,alt_d,azi_d,terrain_alt,camera_polar[0:2,i,j],model_polar[0:2,i,j],blue_camera,blue_model,icloud_camera,icloud_model
        endif

;       Note that when clouds are present the contingency table has a zero value
        icont[1-icloud_camera,1-icloud_model] = icont[1-icloud_camera,1-icloud_model] + 1

        npotl = npotl + 1

        ncloud_camera = ncloud_camera + icloud_camera
        ncloud_model  = ncloud_model  + icloud_model

    endif


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


;   Compare to random analysis having either 50% or the actual percentage of clouds
dirout = lapsdataroot + '/' + 'lapsprd/verif/allsky/stats'
; dirout = '/data/fab/dlaps/projects/roc/hires2/log'
statsfile = 'verif_allsky_anal' + '.' + site + '.' + a9time
openw, 1, dirout + '/' + statsfile

printf,1,'       hits       misses     false_alarms correct_negatives total'
printf,1,        hits,      misses,    false_alarms,correct_negatives,total

printf,1,'     frac_obs     frac_fcst     accuracy   accuracy_random'
printf,1,      frac_obs,    frac_fcst,    accuracy,  accuracy_random

bias = (float(hits) + float(false_alarms)) / (float(hits) + float(misses))

printf,1,FORMAT = '("Bias: ",F7.3,"    ETS: ",F7.3)',bias,ets
printf,1,FORMAT = '(A16,4f10.4," gnuplot")',filetime,frac_obs,frac_fcst,accuracy,accuracy_random
printf,1,FORMAT = '("    solar_alt: ",F7.3)',solar_alt

close,1

exit

end
