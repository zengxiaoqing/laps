
function gen_colormap( )
begin
        ;This procedure explicitly sets the red, green, and blue
        ;intensities for each element in a colormap array.
 
        ; Set the color map
        ; Entry 0 is the background color.
 
        ; Entry 1 is the foreground color.
        ; The continental outlines in the map object use this color
        ; as a default.
        return( (/ (/.000,.000,.000/), \
                (/1.000,1.000,1.000/), \
                (/.000,.000,.000/), \
                (/.020,.020,.020/), \
                (/.040,.040,.040/), \
                (/.060,.060,.060/), \
                (/.080,.080,.080/), \
                (/.100,.100,.100/), \
                (/.120,.120,.120/), \
                (/.140,.140,.140/), \
                (/.160,.160,.160/), \
                (/.180,.180,.180/), \
                (/.200,.200,.200/), \
                (/.220,.220,.220/), \
                (/.240,.240,.240/), \
                (/.260,.260,.260/), \
                (/.280,.280,.280/), \
                (/.300,.300,.300/), \
                (/.320,.320,.320/), \
                (/.340,.340,.340/), \
                (/.360,.360,.360/), \
                (/.380,.380,.380/), \
                (/.400,.400,.400/), \
                (/.420,.420,.420/), \
                (/.440,.440,.440/), \
                (/.460,.460,.460/), \
                (/.480,.480,.480/), \
                (/.500,.500,.500/), \
                (/.520,.520,.520/), \
                (/.540,.540,.540/), \
                (/.560,.560,.560/), \
                (/.580,.580,.580/), \
                (/.600,.600,.600/), \
                (/.620,.620,.620/), \
                (/.640,.640,.640/), \
                (/.660,.660,.660/), \
                (/.680,.680,.680/), \
                (/.700,.700,.700/), \
                (/.720,.720,.720/), \
                (/.740,.740,.740/), \
                (/.760,.760,.760/), \
                (/.780,.780,.780/), \
                (/.800,.800,.800/), \
                (/.820,.820,.820/), \
                (/.840,.840,.840/), \
                (/.860,.860,.860/), \
                (/.880,.880,.880/), \
                (/.900,.900,.900/), \
                (/.920,.920,.920/), \
                (/.940,.940,.940/), \
                (/.960,.960,.960/), \
                (/.980,.980,.980/), \
                (/.990,.990,.990/), \
                (/1.000,1.000,1.000/) /) )
end
;
begin
;
  cdf_file = addfile("static.cdf","r") 
;  cdf_file = addfile("carib.cdf","r") 
;  cdf_file = addfile("datelinetest.cdf","r") 
;  cdf_file = addfile("baja2.cdf","r") ; baja2
;  cdf_file = addfile("taiwan.cdf","r") ; taiwan
;  cdf_file = addfile("alaska.cdf","r") ; alaska
;  cdf_file = addfile("CONUS_10.cdf","r") ; CONUS 10km
;
;;;; Read parameter file as float and as string,
;;;; z1(0) and z1(1) are weights 1 and 2,
;;;; z2(2) and z2(3) are the netCDF variables
;
  z1 = asciiread("albintparms.txt",-1,"float")
  z2 = asciiread("albintparms.txt",-1,"string")
;
;;;; This batch of if statements are necessary since NCL has no
;;;; way to resolve a string variable to use its literal value
;
  if (z2(2) .eq. "a01") then
    var1 = cdf_file->a01(0,0,:,:)
  end if
  if (z2(2) .eq. "a02") then
    var1 = cdf_file->a02(0,0,:,:)
  end if
  if (z2(2) .eq. "a03") then
    var1 = cdf_file->a03(0,0,:,:)
  end if
  if (z2(2) .eq. "a04") then
    var1 = cdf_file->a04(0,0,:,:)
  end if
  if (z2(2) .eq. "a05") then
    var1 = cdf_file->a05(0,0,:,:)
  end if
  if (z2(2) .eq. "a06") then
    var1 = cdf_file->a06(0,0,:,:)
  end if
  if (z2(2) .eq. "a07") then
    var1 = cdf_file->a07(0,0,:,:)
  end if
  if (z2(2) .eq. "a08") then
    var1 = cdf_file->a08(0,0,:,:)
  end if
  if (z2(2) .eq. "a09") then
    var1 = cdf_file->a09(0,0,:,:)
  end if
  if (z2(2) .eq. "a10") then
    var1 = cdf_file->a10(0,0,:,:)
  end if
  if (z2(2) .eq. "a11") then
    var1 = cdf_file->a11(0,0,:,:)
  end if
  if (z2(2) .eq. "a12") then
    var1 = cdf_file->a12(0,0,:,:)
  end if
;
  if (z2(3) .eq. "a01") then
    var2 = cdf_file->a01(0,0,:,:)
  end if
  if (z2(3) .eq. "a02") then
    var2 = cdf_file->a02(0,0,:,:)
  end if
  if (z2(3) .eq. "a03") then
    var2 = cdf_file->a03(0,0,:,:)
  end if
  if (z2(3) .eq. "a04") then
    var2 = cdf_file->a04(0,0,:,:)
  end if
  if (z2(3) .eq. "a05") then
    var2 = cdf_file->a05(0,0,:,:)
  end if
  if (z2(3) .eq. "a06") then
    var2 = cdf_file->a06(0,0,:,:)
  end if
  if (z2(3) .eq. "a07") then
    var2 = cdf_file->a07(0,0,:,:)
  end if
  if (z2(3) .eq. "a08") then
    var2 = cdf_file->a08(0,0,:,:)
  end if
  if (z2(3) .eq. "a09") then
    var2 = cdf_file->a09(0,0,:,:)
  end if
  if (z2(3) .eq. "a10") then
    var2 = cdf_file->a10(0,0,:,:)
  end if
  if (z2(3) .eq. "a11") then
    var2 = cdf_file->a11(0,0,:,:)
  end if
  if (z2(3) .eq. "a12") then
    var2 = cdf_file->a12(0,0,:,:)
  end if
;
  var = var1*z1(0) + var2*z1(1)
;
  lat2d = cdf_file->lat(0,0,:,:)
  lon2d = cdf_file->lon(0,0,:,:)
  varstr = "alb"
  nx = cdf_file->Nx(0)
  xdim = nx - 1
  ny = cdf_file->Ny(0)
  ydim = ny - 1
  rnx = 1.*xdim
  rny = 1.*ydim
  if (rny .ge. rnx) then
    vpheight = .80
    vpy = 1.-((1.-vpheight)/2.)
    vpwidth = (vpheight/rny)*rnx
    vpx = (1.-((vpheight/rny)*rnx))/2.
  else
    vpwidth = .80
    vpx = (1.-vpwidth)/2.
    vpheight = (vpwidth/rnx)*rny
    vpy = 1.-((1.-((vpwidth/rnx)*rny))/2.)
  end if
  latin1 = cdf_file->Latin1
  latin2 = cdf_file->Latin2
  lov = cdf_file->LoV
  proj = cdf_file->grid_type
  projstr = ""
  do n = 0,30
    projstr = projstr + proj(0,n)
  end do
  spac = cdf_file->grid_spacing/1000.
;
; Create an application object.
;
appid = create "wrfsi" appClass defaultapp
    "appUsrDir" : "./"
    "appDefaultParent" : True
end create
;
; Set only one of the following three params to 1, depending on the type
; of output desired (NCGM=metacode file, X11=screen, PS=PostScript file)
;
    NCGM=1
    X11=0
    PS=0
;
    if (NCGM .eq. 1) then
;
; Create an ncgmWorkstation object.
;
        wid = create "wrfsiWork" ncgmWorkstationClass defaultapp
           "wkMetaName" : "./albint.ncgm"
        end create
    else
      if (X11 .eq. 1) then
;
; Create an XWorkstation object.
;
        wid = create "wrfsiWork" xWorkstationClass defaultapp
            "wkPause" : True
        end create
      else
        if (PS .eq. 1) then
;
; Create a PSWorkstation object.
;
          wid = create "wrfsiWork" psWorkstationClass defaultapp
           "wkPSFileName" : "./albint.ps"
          end create
        end if
      end if
    end if
;
; Assign the colormap to the workstation.
setvalues wid
; Generate a colormap.
   "wkColorMap" : gen_colormap()
end setvalues
; Increase workspace memory
setvalues NhlGetWorkspaceObjectId()
   "wsMaximumSize": 67108864
end setvalues
;
difflat = (ny * spac) / 111.1 ; approx latitude degree spacing
difflon = (nx * spac) / 111.1 ; approx longitude degree spacing
gridsp = 1.
if (difflat .gt. 20. .or. difflon .gt. 20.) then
  gridsp = 5.
end if
if (difflat .gt. 40. .or. difflon .gt. 40.) then
  gridsp = 10.
end if
if (projstr .eq. "tangential lambert conformal  " .or.\
    projstr .eq. "secant lambert conformal      ") then
  mapproj = "LambertConformal"
  mapcentlat = 0.0 ; these values are ignored by the Lambert
  mapcentlon = 0.0 ; projection, but need to be set anyway
end if
if (projstr .eq. "mercator                      ") then
  mapproj = "Mercator"
  mapcentlat = 0.0 ; needs to be 0.0 for straight lat lines
  mapcentlon = cdf_file->center_lon
end if
if (projstr .eq. "polar stereographic           ") then
  mapproj = "Stereographic"
  mapcentlat = cdf_file->center_lat ; settings necessary for PS
  mapcentlon = lov                  ; projection
end if
;
mpid = create "mapplot" mapPlotClass wid
;
  "mpProjection" : mapproj
;
  "mpLimitMode" : "Corners"    ; Limit the map view.
  "mpLeftCornerLonF" : lon2d(1,1) 
  "mpLeftCornerLatF" : lat2d(1,1)
  "mpRightCornerLonF" : lon2d(ydim,xdim)
  "mpRightCornerLatF" : lat2d(ydim,xdim)
;
  "mpLambertParallel1F" : latin1
  "mpLambertParallel2F" : latin2
  "mpLambertMeridianF" : lov 
  "mpCenterLonF" : mapcentlon
  "mpCenterLatF" : mapcentlat
  "tfDoNDCOverlay" : True
;
  "mpDataBaseVersion" : "RANGS_GSHHS"
  "mpOutlineBoundarySets" : "GeophysicalAndUSStates"
  "mpGeophysicalLineColor" : "Foreground"
  "mpUSStateLineColor" : "Foreground"
  "mpUSStateLineThicknessF" : 1.25
  "mpGridSpacingF" : gridsp
  "mpGridLineColor" : "Foreground"
  "mpGridLineDashPattern" : 2
  "mpPerimOn" : True
  "mpPerimLineThicknessF" : 1.5
  "vpXF"      : vpx ; Viewport settings
  "vpYF"      : vpy
  "vpWidthF"  : vpwidth
  "vpHeightF" : vpheight
;
end create
;
mpid2 = create "mapplot" mapPlotClass wid
;
; map object strictly to create US state outlines
;
  "mpProjection" : mapproj
;
  "mpLimitMode" : "Corners"    ; Limit the map view.
  "mpLeftCornerLonF" : lon2d(1,1)
  "mpLeftCornerLatF" : lat2d(1,1)
  "mpRightCornerLonF" : lon2d(ydim,xdim)
  "mpRightCornerLatF" : lat2d(ydim,xdim)
;
  "mpLambertParallel1F" : latin1
  "mpLambertParallel2F" : latin2
  "mpLambertMeridianF" : lov
  "mpCenterLonF" : mapcentlon
  "mpCenterLatF" : mapcentlat
  "tfDoNDCOverlay" : True
;
  "mpDataBaseVersion" : "Ncarg4_1"
  "mpOutlineBoundarySets" : "USStates"
  "mpUSStateLineColor" : "Foreground"
  "mpNationalLineColor" : "Foreground"
  "mpGeophysicalLineColor" : "Foreground"
  "mpUSStateLineThicknessF" : 1.25
  "mpOutlineDrawOrder" : "Draw"
  "mpGridSpacingF" : gridsp
  "mpGridLineColor" : "Foreground"
  "mpGridLineDashPattern" : 2
  "mpPerimOn" : True
  "mpPerimLineThicknessF" : 1.5
  "vpXF"      : vpx ; Viewport settings
  "vpYF"      : vpy
  "vpWidthF"  : vpwidth
  "vpHeightF" : vpheight
;
end create
;
; Create a ScalarField object.
;
varfield = create "ScalarField" scalarFieldClass appid
    "sfDataArray" : var
;    "sfXCStartV" : 0 
;    "sfYCStartV" : 0 
;    "sfXCEndV" : xdim 
;    "sfYCEndV" : ydim
    "sfXCStartSubsetV" : 0
    "sfYCStartSubsetV" : 0
    "sfXCEndSubsetV" : xdim - 1
    "sfYCEndSubsetV" : ydim - 1
    "sfMissingValueV" : 1.0E+37
end create
;
if (difflat .lt. 20. .and. difflon .lt. 20.) then
  levsarr = (/0,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,\
              0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40,\
              0.42,0.44,0.46,0.48,0.50,0.52,0.54,0.56,0.58,0.60,\
              0.62,0.64,0.66,0.68,0.70,0.72,0.74,0.76,0.78,0.80,\
              0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98,1.00/)
  colsarr = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,\
              19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,\
              34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,\
              49,50/)
else
  levsarr = (/0.1,0.2,0.3,0.4,0.5,\
              0.6,0.7,0.8,0.9,1.00/)
  colsarr = (/5,10,15,20,25,30,35,40,45,50/)
end if
if (max(var) .lt. 1.) then
  plotmax = 1.
end if
if (max(var) .lt. .9) then
  plotmax = .9
end if
if (max(var) .lt. .8) then
  plotmax = .8
end if
if (max(var) .lt. .7) then
  plotmax = .7
end if
if (max(var) .lt. .6) then
  plotmax = .6
end if
if (max(var) .lt. .5) then
  plotmax = .5
end if
if (max(var) .lt. .4) then
  plotmax = .4
end if
if (max(var) .lt. .3) then
  plotmax = .3
end if
if (max(var) .lt. .2) then
  plotmax = .2
end if
if (max(var) .lt. .1) then
  plotmax = .1
end if
levsp = plotmax / 50.
;
; Create a ContourPlot object.
;
cnid = create "contourplot" contourPlotClass wid
    "cnScalarFieldData":  varfield
    "cnLevelSelectionMode" : "ManualLevels"
    "cnMinLevelValF" : 0. 
    "cnMaxLevelValF" : plotmax
    "cnLevelSpacingF" : levsp
    "cnFillOn" : True
    "cnLinesOn" : False
    "cnLineLabelsOn" : False
    "cnInfoLabelOn" : False
    "pmTickMarkDisplayMode" : "NoCreate"
;    "tmXMajorGrid" : True
;    "tmYMajorGrid" : True
;    "tmXMinorGrid" : True
;    "tmYMinorGrid" : True
    "pmLabelBarDisplayMode" : "ALWAYS"
    "pmLabelBarHeightF" : .05
    "pmLabelBarWidthF" : 0.9
    "pmLabelBarSide" : "bottom"
    "lbOrientation" : "horizontal"
    "lbAutoManage" : True
    "lbLabelFont" : 4
    "lbLabelFontHeightF" : .01
    "lbLabelStride" : 2
    "lbPerimOn" : False
    "lbBoxLinesOn" : 0
    "tiMainFuncCode" : "~"
    "tiMainString" : z2(4)
    "tiMainFont" : 4
    "tiMainFontHeightF" : .015
    "tiMainFontColor" : 1
    "tiMainJust" : "CenterCenter"
    "tiMainOffsetXF" : 0.0
    "tiMainOffsetYF" :	-0.002
    "vpXF"      : vpx ; Viewport settings
    "vpYF"      : vpy
    "vpWidthF"  : vpwidth
    "vpHeightF" : vpheight
;
end create
;
if (difflat .lt. 20. .and. difflon .lt. 20.) then
  boxct = 50
  stride = 2
  filllabels = (/0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,\
              0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40,\
              0.42,0.44,0.46,0.48,0.50,0.52,0.54,0.56,0.58,0.60,\
              0.62,0.64,0.66,0.68,0.70,0.72,0.74,0.76,0.78,0.80,\
              0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98/)
  fillcols = (/3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,\
              19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,\
              34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,\
              49,50,51/)
else
  boxct = 10
  stride = 1
  filllabels = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.00/)
  fillcols = (/5,10,15,20,25,30,35,40,45,50/)
end if
;
lbid = create "LabelBar" labelBarClass wid
;    "lbAutoManage" : False
    "lbAutoManage" : True
    "lbOrientation" : "horizontal"
    "lbLabelFont" : 4
    "lbLabelFontHeightF" : .008
    "lbMonoFillPattern" : True
    "lbBoxLinesOn" : False
    "lbBoxCount" : boxct
    "lbFillColors" : fillcols
    "lbLabelStrings" : filllabels
    "lbLabelStride" : stride
    "lbLabelAlignment" : "InteriorEdges"
    "lbPerimOn" : False
    "vpXF"      : 0.0125
    "vpYF"      : vpy-vpheight
    "vpWidthF"  : 0.975
    "vpHeightF" : 0.04
end create
;
txid = create "TextItems" textItemClass wid
  "txFuncCode" : "~"
  "txString" : "Powered by NCL"
  "txFont" : 4
  "txFontHeightF" : 0.008
  "txFontColor" : 1
end create
;
ianno = NhlAddAnnotation(mpid,txid)
setvalues ianno
;  "amZone":           10
  "amSide":           "bottom"
  "amJust":           "bottomright"
  "amParallelPosF":    0.49
  "amOrthogonalPosF":  0.49
end setvalues
;
; Draw map
;
draw(cnid)
draw(mpid)
draw(mpid2)
frame(wid)
;
end
;
