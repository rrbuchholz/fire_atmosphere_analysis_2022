;=======================================================;
;  plot_regions_overlay_trend.ncl
;=======================================================;
;
; Concepts illustrated:
;   - plot a map of trends with overlayed boxes showing selected regions
; 
;   To run type:
;         ncl plot_regions_overlay_trend.ncl
;
;   Requires:
;         input emission file
;                                  RRB May 5, 2020
;=======================================================;
; load the libraries

begin

;--------------------------------------------
; user input
;--------------------------------------------
  tracer = "CO"
  ; -------------------------------
  ; emission files
  ; -------------------------------
  top_dir = "/MOPITT/V8J/Archive/L3/"
  total_files := systemfunc ("ls "+top_dir+"*08/month/MOP03JM-*.he5")
  monthname = "August"
  years = ispan(2000,2018,1)
  yyyymm = years*100 + 8
  yfrac = yyyymm_to_yyyyfrac(yyyymm,0)

  print(total_files)

  outpath = "/IASI/home/buchholz/"

  ;------------
  ; toggles
  ;------------
  PLOTREGIONS    = True
  if (PLOTREGIONS) then
    mapType2 = "pdf"
    mapName2 = tracer+"_trend_NorthAmerica_"+monthname
  end if ;PLOT
     WITHCO      = True

;--------------------------------------------
; end user input
;--------------------------------------------

;********************************************
; set up
  ; -------------------------------
  ; plotting parameters
  ; -------------------------------
      latmin                         = 20.               ; USA
      latmax                         = 70.
      lonmin                         = -140.
      lonmax                         = -60.

  cr = "~C~"

  ; -------------------------------
  ; Define regions
  ; -------------------------------
   region_names = (/"PNW", "Central USA", "Northeast"/)

   ;                   minlat, maxlat, minlon, maxlon
   region_select = (/(/ 38.,  57.,     -127.,     -110./),\ ;PNW
                     (/ 35.,  49.,     -110.,     -95./),\  ;Central USA
                     (/ 41.,  53.,     -95.,     -74./) /)  ;Northeast


   n_regions = dimsizes(region_select(:,0))

   ; -----------------------------------------------------
   ; collect region outlines for plotting later
     region_lon := region_select(:,0:1)
     region_lat := region_select(:,2:3)


  ; -------------------------------
  ; Load CO and calculate trend
  ; -------------------------------
    fin = addfiles(total_files, "r")
    ListSetType (fin, "join") 
    tracer := fin[:]->RetrievedCOTotalColumnDay_MOP03
      tracer!1 = "longitude"
      tracer!2 = "latitude"
   tracer_array = new((/dimsizes(years), dimsizes(tracer&longitude), dimsizes(tracer&latitude)/), float)

      ; Select years
      if (monthname.eq."August") then
        ; August 2009 is missing
        ; so need to add missing value in array
        ; and remove 2000 and 2001
        tracer_array(2:8,:,:)=tracer(2:8,:,:)
        tracer_array(10:18,:,:)=tracer(9:17,:,:)
      else
        ; Remove 2000 and 2001
        tracer_array(2:20,:,:)=tracer(2:20,:,:)
      end if

      tracer_array!0 = "time"
      tracer_array&time = yyyymm

   mean_co = dim_avg_n_Wrap(tracer_array,0)
      mean_co!0 = "longitude"
      mean_co!1 = "latitude"

   ; calc trend
   rc           = ((regCoef(yfrac,tracer_array(latitude|:,longitude|:,time|:)))/mean_co(latitude|:,longitude|:))*100       
   
   rc@long_name = "regression coefficient (trend)"
   rc@units     = "percent/year"
      rc!0 = "latitude"
      rc&latitude = tracer&latitude
      rc!1 = "longitude"
      rc&longitude = tracer&longitude
 print(rc@tval)

    tval1 = new(dimsizes(rc) , float)
    nxy1  = new(dimsizes(rc) , integer)
    b1    = new(dimsizes(rc) , float)
    unity    = new(dimsizes(rc) , float)
    rc1   = regcoef (yfrac,tracer_array(latitude|:,longitude|:,time|:), tval1,nxy1)
      rc1!0 = "latitude"
      rc1&latitude = tracer&latitude
      rc1!1 = "longitude"
      rc1&longitude = tracer&longitude
    rc1a = (rc1/mean_co(latitude|:,longitude|:))*100


   df   = nxy1-2       ; degrees of freedom: n-2
   b1   = 0.5          ; b must be same dimensions as tval (and df)
   prob = (1 - betainc(df/(df+tval1^2),df/2.0,b1) )
      prob!0 = "latitude"
      prob&latitude = tracer&latitude
      prob!1 = "longitude"
      prob&longitude = tracer&longitude
      prob@_FillValue = -1e+30
    
  ; Select probability limit
   prob_limit = 0.95
 
    rc1b = where(prob.lt.prob_limit, rc1a@_FillValue, rc1a)
      rc1b!0 = "latitude"
      rc1b&latitude = tracer&latitude
      rc1b!1 = "longitude"
      rc1b&longitude = tracer&longitude

   lat_array = conform(rc1b,rc1b&latitude,0)
   lat_array = where(prob.lt.prob_limit, lat_array@_FillValue, lat_array)
   lon_array = conform(rc1b,rc1b&longitude,1)
   lon_array = where(prob.lt.prob_limit, lon_array@_FillValue, lon_array)


printVarSummary(tracer)
printVarSummary(tracer_array)
printVarSummary(rc)
printVarSummary(mean_co)
print(yyyymm)
print(yyyymm+ " "+tracer_array(:,{-150},{40}))
;print(tracer&latitude)

  ;********************************************
  ;-------------------
  ; MAP PLOT
  ;-------------------
  if (PLOTREGIONS) then
    wks2 = gsn_open_wks(mapType2,mapName2)
    ;-------------------
    ; define resources 
    ; to use for plotting
    ;-------------------
    mapres = True
    mapres@tiMainString              = ""                 ; changing main large title
    mapres@gsnLeftString             = ""                 ; change left 
    mapres@gsnRightString             = ""                 ; change left 
    mapres@gsnFrame                  = False              ; do not advance frame
    mapres@gsnDraw		     = False	          ; don't draw it yet
    mapres@mpCenterLonF              =10

    mapres@mpLimitMode               ="LatLon"            ; zoom to Australia
      mapres@mpMinLatF               = latmin
      mapres@mpMaxLatF               = latmax
      mapres@mpMinLonF               = lonmin
      mapres@mpMaxLonF               = lonmax
    ;  mapres@mpCenterLonF            = centerlon

    mapres@vpXF = 0.1
    mapres@vpYF = 0.8
    mapres@vpHeightF            = 0.4          ; set height of plot

    ;-------------------
    ; do the plotting
    ;-------------------

   if (WITHCO) then
     mapres@gsnMaximize           = True             ; make large
     mapres@mpFillOn              = False            ; turn off default background gray
     mapres@mpGeophysicalLineThicknessF = 6
     mapres@tmBorderThicknessF   = 4.0
     mapres@tmYROn = False            ; Turn off right tickmarks.
     mapres@tmXTOn = False            ; Turn off top tickmarks.    

     ; contour levels
     mapres@cnFillOn              = True             ; turn on color
     mapres@cnFillPalette         = "BlueWhiteOrangeRed" ; set color map
     mapres@cnLinesOn             = False            ; turn off contour lines
     mapres@cnLineLabelsOn        = False            ; turn off contour line labels
     mapres@cnFillMode            = "RasterFill"
     mapres@cnLevelSelectionMode  = "ManualLevels"   ; set manual contour levels
     mapres@cnMinLevelValF        =  -1.5           ; set min contour level
     mapres@cnMaxLevelValF        =   1.5           ; set max contour level
     mapres@cnLevelSpacingF       =   0.05            ; set contour interval

     ; label bar
     mapres@lbBoxEndCapStyle    = "TriangleBothEnds"
     mapres@lbBoxLinesOn        = True
     mapres@lbTitleOn           = True
     mapres@lbTitleFontHeightF  = 0.014
     mapres@lbTitleString       = " % per year"
     mapres@lbTitlePosition     = "bottom"
   
     mapres@tiMainString          = ""
     map1 = gsn_csm_contour_map_ce(wks2,rc,mapres)    

    ; Add significance hashes
    pmres = True
      pmres@gsMarkerColor = "gray60"
      pmres@gsMarkerIndex = 16
      pmres@gsMarkerSizeF = 0.004

    dummy = gsn_add_polymarker(wks2, map1, lon_array, lat_array, pmres)


   else 
    map1 = gsn_csm_map(wks2,mapres)                          ; draw a blank map
   end if

    ;-------------------
    ; add names
    ;-------------------
    txres = True
      txres@txFontHeightF = 0.018                     ; font smaller. default big
      txres@txFontThicknessF = 12.
      txres@txFontColor = "black"
      txres@txBackgroundFillColor = "white"
      txres@txPerimSpaceF = 0.2

    do i=0,n_regions-1
      if (region_names(i).eq."Central USA") then 
        dum = gsn_add_text(wks2,map1,region_names(i),(region_select(i,2)+region_select(i,3))/2,\
                           region_select(i,0)-1.5,txres)
      else
        dum = gsn_add_text(wks2,map1,region_names(i),(region_select(i,2)+region_select(i,3))/2+0.8,\
                           region_select(i,0)-1.5,txres)
      end if
    end do

    ; Title not bold
    txres@txFontHeightF = 0.024                     ; font smaller. default big
    gsn_text_ndc(wks2,"a",0.12, 0.87, txres)
    title = monthname+" CO trend"
    gsn_text_ndc(wks2,title,0.52, 0.87, txres)

    draw(map1)

    ;-------------------
    ; draw region boxes
    ;-------------------
    boxres = True
      boxres@gsLineThicknessF = 8

      do j=0,n_regions-1
       boxres@gsLineColor      = "gray20"
       boxres@gsLineDashPattern  = 0 
         regions_lon := (/region_select(j,2), region_select(j,2), \
                          region_select(j,3) , region_select(j,3), \
                          region_select(j,2)/)
         regions_lat := (/region_select(j,0), region_select(j,1), \
                          region_select(j,1) , region_select(j,0), \
                          region_select(j,0)/)
         gsn_polyline(wks2,map1,regions_lon,regions_lat,boxres)
      end do

    ;-------------------
    ; Add significance legend
    ;-------------------
    txres                   = True                     ; polyline mods desired
     txres@txFontHeightF    = 0.018                   ; font smaller. default big
     txres@txJust           = "CenterLeft"
     txres@txPerimOn        = True 
     txres@txFontColor      = "black" 
    pmres = True
    pmres@gsMarkerIndex = 16
    pmres@gsMarkerSizeF = 0.008
    pmres@gsMarkerColor := "gray60"
    pmres@gsFillBackgroundColor = "white"

    gsn_text_ndc(wks2,"    significant trend at ~F33~a ~F21~= 0.05 ",0.118,0.330,txres)
    gsn_polymarker_ndc(wks2,0.125,0.330,pmres)

    frame(wks2)

   end if ; PLOTREGIONS

end

