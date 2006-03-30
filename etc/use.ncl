
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
                (/.500,.500,.500/), \
                (/.500,.500,.500/), \
                (/.980,.910,.570/), \
                (/.940,.880,.550/), \
                (/.890,.830,.510/), \
                (/.800,.740,.460/), \
                (/.600,.560,.350/), \
                (/.500,.520,.290/), \
                (/.600,.825,.600/), \
                (/.500,.775,.500/), \
                (/.400,.725,.400/), \
                (/.300,.675,.300/), \
                (/.200,.625,.200/), \
                (/.100,.575,.100/), \
                (/.000,.500,.000/), \
                (/.000,.400,.000/), \
                (/.160,.390,.660/), \
                (/.400,.200,.400/), \
                (/.500,.300,.500/), \
                (/.600,.400,.600/), \
                (/.700,.500,.700/), \
                (/.500,.300,.300/), \
                (/.750,.500,.500/), \
                (/.850,.700,.700/), \
                (/.950,.950,.950/), \
                (/.000,.500,.000/) /) )
end
;
begin

  cdf_file = addfile("static.cdf","r") 
;  cdf_file = addfile("baja2.cdf","r") ; baja2
;  cdf_file = addfile("carib.cdf","r") ; caribbean
;  cdf_file = addfile("taiwan.cdf","r") ; taiwan
;  cdf_file = addfile("alaska.cdf","r") ; alaska
;  cdf_file = addfile("capecod.cdf","r") ; cape cod
;
  var = cdf_file->use(0,0,:,:)
  lat2d = cdf_file->lat(0,0,:,:)
  lon2d = cdf_file->lon(0,0,:,:)
  varstr = "use"
  nx = cdf_file->Nx(0)
  xdim = nx - 1
  ny = cdf_file->Ny(0)
  ydim = ny - 1
  rnx = 1.*xdim
  rny = 1.*ydim
  if (rny .ge. rnx) then
    vpheight = .75
    vpy = 1.-((1.-vpheight)/2.)
    vpwidth = (vpheight/rny)*rnx
    vpx = (1.-((vpheight/rny)*rnx))/2.
  else
    vpwidth = .75
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
           "wkMetaName" : "./use.ncgm"
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
           "wkPSFileName" : "./use.ps"
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
   "wsMaximumSize": 33554432
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
;  "mpLeftCornerLonF" : lon2d(0,0) 
;  "mpLeftCornerLatF" : lat2d(0,0) 
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
  "mpGeophysicalLineColor" : "Background"
  "mpUSStateLineColor" : "Background"
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
;  "mpLeftCornerLonF" : lon2d(0,0)
;  "mpLeftCornerLatF" : lat2d(0,0)
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
  "mpUSStateLineColor" : "Background"
  "mpNationalLineColor" : "Background"
  "mpGeophysicalLineColor" : "Background"
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
fillindices = ispan(3,26,1)
minlv = 1
maxlv = 24
levsp = 1
;
; Create a ContourPlot object.
;
cnid = create "contourplot" contourPlotClass wid
    "cnScalarFieldData":  varfield
    "cnFillColors": fillindices
    "cnLevelSelectionMode" : "ManualLevels"
    "cnMinLevelValF" : minlv 
    "cnMaxLevelValF" : maxlv
    "cnLevelSpacingF" : levsp 
    "cnRasterModeOn" : True
    "cnFillOn" : True
    "cnLinesOn" : False
    "cnLineLabelsOn" : False
    "cnInfoLabelOn" : False
    "pmTickMarkDisplayMode" : "NoCreate"
;    "tmXMajorGrid" : True
;    "tmYMajorGrid" : True
;    "tmXMinorGrid" : True
;    "tmYMinorGrid" : True
;    "tmXTMinorPerMajor" : 1
;    "tmXBMinorPerMajor" : 1
;    "tmYLMinorPerMajor" : 1
;    "tmYRMinorPerMajor" : 1
    "tiMainString" : var@long_name
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
lgid = create "Legend" legendClass wid
    "lgAutoManage" : False
    "lgOrientation" : "Vertical"
    "lgItemCount" : 16
    "lgLineColors" : (/4,5,6,7,8,9,10,11,12,13,14,15,16,3,17,18/)
    "lgLabelStrings" : (/"Sand","Loamy Sand","Sandy Loam","Silt Loam",\
        "Silt","Loam","Sandy Clay Loam","Silty Clay Loam","Clay Loam",\
	"Sandy Clay","Silty Clay","Clay","Org. Materials","Water",\
        "Bedrock","Other (land-ice)"/)
    "lgLabelFont" : 4
    "lgLabelFontHeightF" : 0.005
    "lgLabelOffsetF" : 0.03
    "lgLineThicknessF" : 4.0
    "lgMonoDashIndex" : True
    "vpXF"	: vpx + vpwidth
    "vpYF"	: vpy
    "vpWidthF"	: 0.03
    "vpHeightF" : vpheight
end create
;
lbid = create "LabelBar" labelBarClass wid
    "lbAutoManage" : False
    "lbLabelFont" : 4
    "lbLabelFontHeightF" : .008
    "lbMonoFillPattern" : True
    "lbBoxCount" : 9
    "lbBoxLinesOn" : True
    "lbBoxCount" : 24
    "lbFillColors" : (/3,4,5,6,7,8,9,10,11,12,13,14,15,\
        16,17,18,19,20,21,22,23,24,25,26,27/)
    "lbLabelAlignment" : "BoxCenters"
    "lbLabelStrings" : (/"Urban","Dry Cropland","Irr Cropland",\
        "Mixed Dry/Irr","Crop/Grassland","Crop/Woodland","Grassland",\
        "Shrubland","Mixed Shr/Grs","Savanna","D. Broad Forest",\
        "D. Needle Forest","E. Broad Forest","E. Needle Forest",\
        "Mixed Forest","Water","Herb Wetland","Wood Wetland",\
        "Barren","Herb Tundra","Wood Tundra","Mixed Tundra",\
        "Bare Gnd Tundra","Snow or Ice"/)
    "vpXF"      : vpx + vpwidth
    "vpYF"      : vpy 
    "vpWidthF"  : 0.04
    "vpHeightF" : vpheight
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
draw(lbid)
frame(wid)
;
end
;
