
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
                (/.200,.200,.200/), \
                (/.225,.225,.225/), \
                (/.250,.250,.250/), \
                (/.275,.275,.275/), \
                (/.300,.300,.300/), \
                (/.325,.325,.325/), \
                (/.350,.350,.350/), \
                (/.375,.375,.375/), \
                (/.400,.400,.400/), \
                (/.425,.425,.425/), \
                (/.450,.450,.450/), \
                (/.475,.475,.475/), \
                (/.500,.500,.500/), \
                (/.525,.525,.525/), \
                (/.550,.550,.550/), \
                (/.575,.575,.575/), \
                (/.600,.600,.600/), \
                (/.160,.390,.660/), \
                (/.090,.530,.230/), \
                (/.109,.535,.235/), \
                (/.127,.540,.240/), \
                (/.145,.545,.245/), \
                (/.163,.550,.250/), \
                (/.182,.555,.255/), \
                (/.200,.560,.260/), \
                (/.218,.565,.265/), \
                (/.236,.570,.270/), \
                (/.255,.575,.275/), \
                (/.273,.580,.280/), \
                (/.291,.585,.285/), \
                (/.309,.590,.290/), \
                (/.327,.595,.295/), \
                (/.346,.600,.300/), \
                (/.364,.605,.305/), \
                (/.382,.610,.310/), \
                (/.400,.615,.315/), \
                (/.419,.620,.320/), \
                (/.437,.625,.325/), \
                (/.455,.630,.330/), \
                (/.474,.635,.335/), \
                (/.492,.640,.340/), \
                (/.510,.645,.345/), \
                (/.528,.650,.350/), \
                (/.546,.655,.355/), \
                (/.565,.660,.360/), \
                (/.583,.665,.365/), \
                (/.601,.670,.370/), \
                (/.618,.675,.375/), \
                (/.636,.680,.380/), \
                (/.653,.685,.385/), \
                (/.670,.690,.390/), \
                (/.665,.675,.385/), \
                (/.660,.660,.380/), \
                (/.655,.645,.375/), \
                (/.650,.630,.370/), \
                (/.645,.615,.365/), \
                (/.640,.600,.360/), \
                (/.635,.585,.355/), \
                (/.630,.570,.350/), \
                (/.625,.555,.345/), \
                (/.620,.540,.340/), \
                (/.615,.525,.335/), \
                (/.610,.510,.330/), \
                (/.605,.495,.325/), \
                (/.600,.480,.320/), \
                (/.595,.465,.315/), \
                (/.590,.450,.310/), \
                (/.585,.435,.305/), \
                (/.580,.420,.300/), \
                (/.575,.405,.295/), \
                (/.570,.390,.290/), \
                (/.565,.375,.285/), \
                (/.560,.360,.280/), \
                (/.555,.345,.275/), \
                (/.550,.330,.270/), \
                (/.545,.315,.265/), \
                (/.540,.300,.260/), \
                (/.535,.285,.255/), \
                (/.530,.270,.250/), \
                (/.522,.260,.243/), \
                (/.515,.250,.235/), \
                (/.507,.240,.227/), \
                (/.500,.230,.220/), \
                (/.506,.235,.231/), \
                (/.512,.240,.243/), \
                (/.518,.244,.254/), \
                (/.523,.249,.265/), \
                (/.529,.254,.276/), \
                (/.535,.259,.288/), \
                (/.541,.264,.299/), \
                (/.546,.268,.310/), \
                (/.552,.273,.322/), \
                (/.558,.277,.333/), \
                (/.564,.282,.344/), \
                (/.569,.286,.355/), \
                (/.575,.291,.366/), \
                (/.581,.296,.378/), \
                (/.586,.301,.389/), \
                (/.592,.305,.400/), \
                (/.598,.310,.412/), \
                (/.604,.315,.423/), \
                (/.609,.319,.434/), \
                (/.615,.324,.445/), \
                (/.621,.329,.456/), \
                (/.627,.333,.468/), \
                (/.632,.338,.479/), \
                (/.638,.342,.490/), \
                (/.644,.347,.502/), \
                (/.650,.352,.513/), \
                (/.655,.356,.524/), \
                (/.661,.361,.535/), \
                (/.669,.366,.547/), \
                (/.676,.371,.558/), \
                (/.683,.375,.569/), \
                (/.690,.380,.580/), \
                (/.688,.382,.582/), \
                (/.686,.383,.584/), \
                (/.683,.384,.586/), \
                (/.681,.386,.588/), \
                (/.679,.388,.590/), \
                (/.677,.390,.592/), \
                (/.675,.391,.594/), \
                (/.673,.393,.595/), \
                (/.671,.395,.597/), \
                (/.669,.396,.599/), \
                (/.666,.397,.601/), \
                (/.664,.399,.603/), \
                (/.662,.400,.605/), \
                (/.660,.402,.607/), \
                (/.658,.404,.609/), \
                (/.655,.405,.610/), \
                (/.653,.406,.612/), \
                (/.651,.408,.614/), \
                (/.649,.410,.616/), \
                (/.646,.411,.618/), \
                (/.644,.413,.620/), \
                (/.642,.415,.622/), \
                (/.640,.417,.623/), \
                (/.638,.418,.625/), \
                (/.636,.419,.627/), \
                (/.634,.421,.629/), \
                (/.632,.423,.631/), \
                (/.629,.424,.633/), \
                (/.627,.425,.635/), \
                (/.625,.427,.637/), \
                (/.623,.428,.638/), \
                (/.620,.430,.640/), \
                (/.629,.440,.650/), \
                (/.638,.450,.660/), \
                (/.647,.460,.669/), \
                (/.656,.470,.679/), \
                (/.666,.480,.689/), \
                (/.675,.490,.699/), \
                (/.684,.500,.709/), \
                (/.693,.510,.718/), \
                (/.702,.520,.727/), \
                (/.711,.530,.737/), \
                (/.720,.540,.747/), \
                (/.729,.550,.756/), \
                (/.738,.560,.766/), \
                (/.747,.570,.776/), \
                (/.756,.580,.786/), \
                (/.765,.590,.795/), \
                (/.774,.600,.805/), \
                (/.783,.610,.815/), \
                (/.792,.620,.824/), \
                (/.801,.630,.834/), \
                (/.811,.640,.844/), \
                (/.820,.650,.854/), \
                (/.829,.660,.863/), \
                (/.838,.670,.873/), \
                (/.847,.680,.883/), \
                (/.856,.690,.892/), \
                (/.865,.700,.901/), \
                (/.874,.710,.911/), \
                (/.883,.720,.921/), \
                (/.892,.730,.931/), \
                (/.901,.740,.941/), \
                (/.910,.750,.950/), \
                (/.913,.757,.947/), \
                (/.916,.765,.943/), \
                (/.919,.772,.939/), \
                (/.921,.779,.935/), \
                (/.924,.786,.931/), \
                (/.927,.794,.928/), \
                (/.930,.801,.924/), \
                (/.933,.808,.920/), \
                (/.936,.815,.916/), \
                (/.939,.822,.913/), \
                (/.942,.829,.909/), \
                (/.944,.836,.905/), \
                (/.947,.844,.901/), \
                (/.949,.851,.898/), \
                (/.952,.858,.894/), \
                (/.955,.865,.890/), \
                (/.957,.873,.887/), \
                (/.959,.880,.883/), \
                (/.962,.887,.879/), \
                (/.966,.894,.875/), \
                (/.969,.901,.872/), \
                (/.972,.909,.868/), \
                (/.975,.916,.864/), \
                (/.978,.923,.860/), \
                (/.981,.930,.857/), \
                (/.984,.937,.853/), \
                (/.987,.944,.849/), \
                (/.989,.951,.845/), \
                (/.992,.959,.842/), \
                (/1.000,0.980,0.830/), \
                (/1.000,0.983,0.857/), \
                (/1.000,0.985,0.885/), \
                (/1.000,0.988,0.913/), \
                (/1.000,0.990,0.940/), \
                (/1.000,0.993,0.955/), \
                (/1.000,0.995,0.970/), \
                (/1.000,0.996,0.980/), \
                (/1.000,0.997,0.982/), \
                (/1.000,0.998,0.984/), \
                (/1.000,0.998,0.986/), \
                (/1.000,0.998,0.988/), \
                (/1.000,0.998,0.990/), \
                (/1.000,0.999,0.992/), \
                (/1.000,0.999,0.994/), \
                (/1.000,0.999,0.996/), \
                (/1.000,0.999,0.998/), \
                (/1.000,1.000,1.000/) /) )
end
;
begin

  cdf_file = addfile("static.cdf","r") 
;  cdf_file = addfile("carib.cdf","r") ; caribbean
;  cdf_file = addfile("baja2.cdf","r") ; baja2
;  cdf_file = addfile("taiwan.cdf","r") ; taiwan
;  cdf_file = addfile("alaska.cdf","r") ; alaska
;  cdf_file = addfile("CONUS_10.cdf","r") ; CONUS_10
;  cdf_file = addfile("capecod.cdf","r") ; cape cod
;
  var = cdf_file->avg(0,0,:,:)
  lat2d = cdf_file->lat(0,0,:,:)
  lon2d = cdf_file->lon(0,0,:,:)
  varstr = "avg"
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
           "wkMetaName" : "./avg.ncgm"
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
           "wkPSFileName" : "./avg.ps"
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
  "mpFillAreaSpecifiers" : (/"USStatesLand"/)
  "mpSpecifiedFillColors" : (/0/)
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
if (min(var) .eq. max(var)) then
  levsarr = min(var)
  colsarr = 1
  levsp = 1
  lbstr = 1
  paraladj = 0.5
else
  nlevs = 50
  npts = (nx * 1.0) * (ny * 1.0)
  diff = floattoint(max(var)-min(var))
  levsp = floattoint(((max(var)-min(var))/nlevs)/10.) * 10
  firstlev = .01
  if (levsp .eq. 0) then
    levsp = 1
  end if
  colstep = floattoint(230/nlevs) - 1
  minlv = floattoint(min(var))
  maxlv = floattoint(max(var))
  startpt = 19 ; 19 is the colormap value for water
  startlv = 0
  if (minlv.lt.0) then
    do while (startlv.gt.minlv)
      startlv = startlv - levsp
      startpt = startpt - 1
    end do
  end if
  if (startpt.lt.0) then
    print("COLORMAP ERROR - startpt less than zero")
  end if
  testlevsarr = ispan(minlv,maxlv,levsp) * 1.0
  colsarr = ispan(startpt,230,colstep)
  colsarr(0) = 0
  colsarr(1) = 19
  do while(dimsizes(colsarr) .le. dimsizes(testlevsarr))
    delete(testlevsarr)
    levsp = levsp + 1
    testlevsarr = ispan(minlv,maxlv,levsp) * 1.0
  end do
  lbstr = 5
  if (dimsizes(colsarr) .le. 15) then
    lbstr = 2
  end if
  paraladj = 0.5 - ((vpwidth/dimsizes(colsarr)) / 2.)
  levsarr = testlevsarr
  levsarr(1) = firstlev
end if
;
; Create a ContourPlot object.
;
cnid = create "contourplot" contourPlotClass wid
    "cnScalarFieldData":  varfield
    "cnLevelSelectionMode" : "ExplicitLevels"
    "cnLevels" : levsarr
    "cnFillColors" : colsarr
    "cnFillOn" : True
;    "cnRasterModeOn" : True
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
    "pmLabelBarWidthF" : vpwidth
    "pmLabelBarSide" : "bottom"
    "pmLabelBarOrthogonalPosF" : 0.0
    "pmLabelBarParallelPosF" : paraladj
    "lbOrientation" : "horizontal"
    "lbAutoManage" : False
    "lbLabelFont" : 4
    "lbLabelFontHeightF" : .01
    "lbLabelStride" : lbstr
    "lbTitleFont" : 22  ; Helvetica-bold
    "lbTitleFontColor" : 1
    "lbPerimOn" : False
    "lbBoxLinesOn" : 0
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
