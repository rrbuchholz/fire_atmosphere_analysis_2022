;  modis_seasonal_cycles.ncl 
;=======================================================;
;
; Concepts illustrated:
;   - For use with TERRA/MODIS extracted seasonal cycles
;
;   Files needed:
;         -MODIS extracted seasonal cycles
;
;   To run type:
;         > ncl modis_seasonal_cycles.ncl
;                                  RRB Novemner 18, 2019
;=======================================================;
; load the libraries
load "./calc_seasonal_cycle.ncl"

begin

;--------------------------------------------
; user input
;--------------------------------------------
location = "Northeast"
plot_title = "Northeast"

  ; -------------------------------
  ; MODIS
  ; -------------------------------
  timeseries_files = systemfunc ("ls /IASI/home/buchholz/MODIS/month_files/modis_MOD08_M3L3_3d.nc")

print(timeseries_files)


  ;------------
  ; toggles
  ;------------
  PLOT           = True
    plttype      = "x11"
    miny = 0.05
    maxy = 0.50
  
;--------------------------------------------
; end user input
;--------------------------------------------
;
;--------------------------------------------
; set up
;--------------------------------------------
  ; -------------------------------
  ; REGION PARAMETERS
  ; -------------------------------
   region_names = (/"PNW", "CentralUSA", "Northeast"/)

   ;                   minlat, maxlat, minlon, maxlon
   region_select = (/(/ 38.,  57.,     -127.,     -110./),\ ;PNWb
                     (/ 35.,  49.,     -110.,     -95./),\  ;Central Mountain USA
                     (/ 41.,  53.,     -95.,     -74./) /)  ;Northeast

   ind_region = ind(region_names.eq.location)

  ; -------------------------------
  ; OUTPUT RELATED SETUP
  ; -------------------------------
  if (PLOT) then
    pltdir       = "./"
    pltname      = "modis_seasonal_cycle" + location
  end if


;--------------------------------------------
; load file and extract
;--------------------------------------------
  ; -------------------------------
  ; MODIS
  ; -------------------------------
     print("Collecting MODIS")
     fin = addfile(timeseries_files, "r")            ; read in
     timeseries_read  = fin->DTDB_CMB_MEAN
        timeseries_read!0 = "time"
        timeseries_read!1 = "lat"
        timeseries_read!2 = "lon"
        timeseries_read@_FillValue = -9999

  
     ; landsea mask
     a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
      lsdata = a->LSMASK
      ;lsm  = landsea_mask(lsdata,timeseries_read&lat,timeseries_read&lon)
      ;timeseries_read = mask(timeseries_read,lsm.eq.0,False)

     timeseries_region = timeseries_read(:,{region_select(ind_region,0):region_select(ind_region,1)},\
                                            {region_select(ind_region,2):region_select(ind_region,3)})

     timeseries_dummy = dim_avg_n(dim_avg_n(timeseries_region,1),1)
     time_dummy = timeseries_region&time
     time_dummy@units = "days since -4712-01-01"
     date_dummy = cd_calendar(time_dummy,1)

     ind_range= ind(date_dummy.ge.200200.and.date_dummy.le.201812)
     timeseries = timeseries_dummy(ind_range)
     time = time_dummy(ind_range)

     time@units = "days since -4712-01-01"
     date = cd_calendar(time,-1)
     modis_yfrac = yyyymm_to_yyyyfrac(date, 0.5)

  ; -------------------------------
  ; PROCESS
  ; -------------------------------

  ; -------------------------------
  ; seasonal_cycles
  ; -------------------------------
   timeseries!0 = "time"
   timeseries&time = date

  modis_seas_avg_a = calc_seasonal_cycle(timeseries,date,2002,2011, 12)
  modis_seas_avg_b = calc_seasonal_cycle(timeseries,date,2012,2018, 12)

  print("2002-2011 AVG --- STD DEV")
  print(modis_seas_avg_a(:,0) + " --- " + modis_seas_avg_a(:,1))
  print("2012-2018 AVG --- STD DEV")
  print(modis_seas_avg_b(:,0) + " --- " + modis_seas_avg_b(:,1))


  ; -------------------------------
  ; CALCULATE t-test SIGNIFICANCE
  ; -------------------------------
  var_a = modis_seas_avg_a(:,1)*modis_seas_avg_a(:,1)
  var_b = modis_seas_avg_b(:,1)*modis_seas_avg_b(:,1)
  num_a = modis_seas_avg_a(:,0)*0
  num_a = 10
  num_b = modis_seas_avg_b(:,0)*0
  num_b = 7

  df = num_a+num_b-2
  pooled_sd = sqrt((var_a+var_b)/2)
  mean1 = modis_seas_avg_a(:,0)
  mean2 = modis_seas_avg_b(:,0)
  cohen_d = abs(mean1-mean2)/pooled_sd

  t_test_result = ttest(modis_seas_avg_a(:,0),var_a,num_a,modis_seas_avg_b(:,0),var_b,num_b,False,True)
  print(t_test_result)

  ;---- 95 % t= 2.447, 90% t = 1.934
  CI_lb = ((modis_seas_avg_a(:,0)-modis_seas_avg_b(:,0))-2.447*sqrt(var_a/num_a+var_b/num_b))
  CI_ub = ((modis_seas_avg_a(:,0)-modis_seas_avg_b(:,0))+2.447*sqrt(var_a/num_a+var_b/num_b))
  mean_diff = (modis_seas_avg_a(:,0)-modis_seas_avg_b(:,0))

  print("*****************************")
  print("August t-test results:")
  print("p-value =   " + sprintf("%9.3f",t_test_result(0,7)) )
  print("t-value =   " + sprintf("%9.3f", t_test_result(1,7)) )
  print("df =        " + df(7))
  print("Cohen's d = " + sprintf("%9.3f", cohen_d(7)) )
  print("95%CI =     " + sprintf("%9.3f", CI_lb(7))+ ";" +sprintf("%9.3f", CI_ub(7)) )
  print("Mean diff = " + sprintf("%9.3f", mean_diff(7)) )
  print("*****************************")

sig_a = where(t_test_result(0,:).lt.(0.05), modis_seas_avg_a(:,0), modis_seas_avg_a@_FillValue)
sig_b = where(t_test_result(0,:).lt.(0.05), modis_seas_avg_b(:,0), modis_seas_avg_b@_FillValue)
sig_a01 = where(t_test_result(0,:).lt.(0.01), modis_seas_avg_a(:,0), modis_seas_avg_a@_FillValue)
sig_b01 = where(t_test_result(0,:).lt.(0.01), modis_seas_avg_b(:,0), modis_seas_avg_b@_FillValue)

;********************************************
;--------------------------------------------
; plot
;--------------------------------------------
if (PLOT) then

  wks = gsn_open_wks(plttype,pltname)        
  gsn_define_colormap(wks,"BlAqGrYeOrRe")  
   
  res                       = True                      ; plot mods desired
    res@gsnDraw		    = False	                ; don't draw it yet
    res@gsnFrame            = False                     ; don't advance frame
    res@gsnMaximize         = True

      res@trYMaxF             = maxy
      res@trYMinF             = miny
    res@tmYLMode        = "Manual"	
      res@tmYLTickStartF  = 0
      res@tmYLTickEndF    = maxy
      res@tmYLTickSpacingF= 0.2

    res@trXMinF             = -0.5
    res@trXMaxF             = 11.5
    res@tmXBMode              = "Explicit"
      res@tmXBValues           = modis_seas_avg_a&month(1::2)
      res@tmXBLabels           = (/"Feb", "Apr", "Jun",\
                                   "Aug", "Oct", "Dec"/)
      res@tmXBLabelFontHeightF = 0.048
      res@tmYLLabelFontHeightF = 0.048
      res@tiYAxisFontHeightF  = 0.048
      res@tiXAxisFontHeightF  = 0.048

    res@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
    res@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticks inward
    res@tmXBMajorOutwardLengthF = 0.0               
    res@tmXBMinorOutwardLengthF = 0.0        

    res@vpWidthF             = 1.0
    res@vpHeightF            = 0.5          ; set height of plot

    res@tiMainString        = ""
    res@tiXAxisString       = "Month"
    res@tiYAxisString       = "AOD at 550 nm"
    res@gsnLeftString       = plot_title

    res@pmLegendDisplayMode      = "Always"        ; turn on legend
    res@lgPerimOn                = False           ; Turn off perimeter
    res@pmLegendWidthF           = 0.15            ; Change width and
    res@pmLegendHeightF          = 0.20            ; height of legend
    res@lgLabelFontHeightF       = 0.042
    res@lgJustification          = "Centerleft"
    res@pmLegendParallelPosF     = 0.60

  ; -------------------------------
  ; timeseries - weekly averages and climatology
  ; -------------------------------
    res@xyExplicitLegendLabels   = "2002-2011"
    res@pmLegendOrthogonalPosF   = -1.42
    res@xyLineColor         = "black"
    res@xyLineThicknessF    = 10.5
    res@xyDashPatterns      =2
  plot1 = gsn_csm_xy(wks,modis_seas_avg_a&month,modis_seas_avg_a(:,0),res)


    res@xyExplicitLegendLabels   = "2012-2018"
    res@pmLegendOrthogonalPosF   = -1.53
    res@gsnLeftString       = ""
    res@xyLineColor         = "red3"
    res@xyDashPatterns      =3
  plot2 = gsn_csm_xy(wks,modis_seas_avg_b&month,modis_seas_avg_b(:,0),res)
  overlay(plot1, plot2)

  res@pmLegendDisplayMode      = "Never"        ; turn on legend

  ; -------------------------------
  ; add std dev
  ; -------------------------------

  ; -------------------------------
  ; Define a polygon for std-dev
  ; -------------------------------
   month_time = modis_seas_avg_a&month
   ya1    = new( (/2*dimsizes(month_time)/), float )
   yb1    = new( (/2*dimsizes(month_time)/), float )
   do k=0,dimsizes(month_time)-1
     ya1(k) = month_time(k)
     yb1(k) = tofloat((modis_seas_avg_a(k,0)+modis_seas_avg_a(k,1)))
     yb1(2*dimsizes(month_time)-1-k) = tofloat((modis_seas_avg_a(k,0)-modis_seas_avg_a(k,1)))
     ya1(2*dimsizes(month_time)-1-k) = month_time(k)
   end do
    gsres                          = True                ; poly res
    gsres@gsFillColor            = "gray60"                ; color chosen
    gsres@gsFillOpacityF         = 0.3                 
    dum1a = gsn_add_polygon(wks,plot1,ya1,yb1,gsres)

  ; -------------------------------
  ; Define a polygon for std-dev
  ; -------------------------------
   ya2    = new( (/2*dimsizes(month_time)/), float )
   yb2    = new( (/2*dimsizes(month_time)/), float )
   do k=0,dimsizes(month_time)-1
     ya2(k) = month_time(k)
     yb2(k) = tofloat((modis_seas_avg_b(k,0)+modis_seas_avg_b(k,1)))
     yb2(2*dimsizes(month_time)-1-k) = tofloat((modis_seas_avg_b(k,0)-modis_seas_avg_b(k,1)))
     ya2(2*dimsizes(month_time)-1-k) = month_time(k)
   end do
    gsres                          = True                ; poly res
    gsres@gsFillColor            = "indianred2"                ; color chosen
    gsres@gsFillOpacityF         = 0.3                 
    dum2a = gsn_add_polygon(wks,plot2,ya2,yb2,gsres)

  ; -------------------------------
  ; plot significant points
  ; -------------------------------
    res@xyMarkLineModes = "Markers"
    res@xyMarkerColors = "black"
    res@xyMarkers = 16
    res@xyMarkerSizeF    = 25.5
  plot1s = gsn_csm_xy(wks,modis_seas_avg_a&month,sig_a,res)
    res@xyMarkerSizeF    = 25.5
  ;plot1s01 = gsn_csm_xy(wks,modis_seas_avg_a&month,sig_a01,res)
  overlay(plot1, plot1s)
  ;overlay(plot1, plot1s01)

    res@xyMarkerColors = "red3"
    res@xyMarkerSizeF    = 25.5
  plot2s = gsn_csm_xy(wks,modis_seas_avg_b&month,sig_b,res)
    res@xyMarkerSizeF    = 25.5
  ;plot2s01 = gsn_csm_xy(wks,modis_seas_avg_a&month,sig_b01,res)
  overlay(plot1, plot2s)
  ;overlay(plot1, plot2s01)

  ; add label
    txres                   = True                     ; polyline mods desired
     txres@txFontHeightF    = 0.032                   ; font smaller. default big
     txres@txJust           = "CenterCenter"
     txres@txPerimOn        = False
     txres@txFontColor      = "black" 
    pmres = True
    pmres@gsMarkerIndex = 16
    pmres@gsMarkerSizeF = 0.012
    pmres@gsMarkerColor := "black"
    gsn_polymarker_ndc (wks,0.72,0.60,pmres)
    gsn_text_ndc(wks,"p < 0.05",0.816,0.60,txres)
    ;pm01res = True
    ;pm01res@gsMarkerIndex = 16
    ;pm01res@gsMarkerSizeF = 0.020
    ;pm01res@gsMarkerColor := "black"
    ;gsn_polymarker_ndc (wks,0.15,0.69,pm01res)
    ;gsn_text_ndc(wks,"p < 0.01",0.214,0.69,txres)


  ; -------------------------------
  ; ADD CUSTOM TEXT
  ; -------------------------------
  ;drawNDCGrid(wks)
  ; STATIONS
  txres                   = True                ; polyline mods desired
  txres@txFontHeightF = 0.022                   ; font smaller. default big
  txres@txJust = "CenterLeft"                   ; font smaller. default big
  txres@txFontColor = "black"
  ;gsn_text_ndc(wks,"MODIS 2002-2011",0.68,0.69,txres)
  txres@txFontColor = "red3"
  ;gsn_text_ndc(wks,"MODIS 2012-2018",0.68,0.66,txres)

  ; -------------------------------
  ; DRAW IT OUT
  ; -------------------------------
  draw(plot1)
  frame(wks)

end if



end

