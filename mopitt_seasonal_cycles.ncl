;=======================================================;
;  mopitt_seasonal_cycles.ncl 
;=======================================================;
;
; Concepts illustrated:
;   - For use with MOPITT CO version 8, region avg data
;   - loading CO data
;   - monthly mean values
;   - whole dataset average month values
;
;   Files needed:
;         -MOPITT subset of level 2, daily hdf files
;
;   To run type:
;         > ncl mopitt_seasonal_cycles.ncl
;                                  RRB Octoeber 20, 2019
;=======================================================;
; load the libraries
load "./calc_seasonal_cycle.ncl"

begin

;--------------------------------------------
; user input
;--------------------------------------------
location = "PNW"
plot_title = "PNW"
obstype = "column"
y_end = 2018

  ; -------------------------------
  ; MOPITT
  ; -------------------------------
  timeseries_files = systemfunc ("ls /IASI/home/buchholz/MOPITT_subset/regions_v8/"+location+"_*monthavg_VMR.nc")
print(timeseries_files)

  ;------------
  ; toggles
  ;------------
  PLOT           = True
    plttype      = "x11"
    if (obstype.eq."surface") then
      y_title = "Surface CO (ppb)"
      miny=110
      maxy=180
    else if (obstype.eq."column")
      y_title = "Column CO (ppb)"
      miny=80
      maxy=130
    end if
    end if

  
;--------------------------------------------
; end user input
;--------------------------------------------
;
;--------------------------------------------
; set up
;--------------------------------------------
  ; -------------------------------
  ; OUTPUT RELATED SETUP
  ; -------------------------------
  if (PLOT) then
    pltdir       = "./"
    pltname      = "mopitt_seasonal_cycle_" +obstype+"_"+ location
  end if

;--------------------------------------------
; load file and extract
;--------------------------------------------
  ; -------------------------------
  ; MOPITT
  ; -------------------------------
  ; read data from region
  fin = addfiles(timeseries_files, "r")      ; read in
    ListSetType (fin, "cat")                 ; concatenate

    if (obstype.eq."surface") then
      ; Surface vmr
      ;timeseries_dummy  = fin[:]->AvgRetrievedCOMixingRatioProfile
      timeseries_dummy  = fin[:]->AvgSurfaceVMR
      global_trend = -0.5
      scale = 1
    else if (obstype.eq."column")
      ; Column vmr
      timeseries_dummy  = fin[:]->RetrievedX_CORegionStats
      scale = 1e-09
      global_trend = -0.5*1e-09
    end if
    end if

    dfs_dummy = fin[:]->AvgDegreesofFreedomforSignal
    ;printVarSummary(timeseries_dummy)

    time_dummy        = fin[:]->time
     date_dummy  = str_split_by_length(tostring(time_dummy),4)
     years_dummy        = tointeger(date_dummy(:,0))
     ind_start = ind(years_dummy.gt.2001)
     ind_end = ind(years_dummy.le.y_end)


    time       = time_dummy(ind_start(0):ind_end(dimsizes(ind_end)-1))
     mopitt_date  = str_split_by_length(tostring(time),4)
     years        = tointeger(mopitt_date(:,0))
     months       = tointeger(mopitt_date(:,1))
     mopitt_yfrac = yyyymm_to_yyyyfrac(time,0.5)
    ; For temp surface layers
    ;timeseries = dim_avg_n(timeseries_dummy(ind_start(0):ind_end(dimsizes(ind_end)-1),0:1),1)
    ; For column vals and true surface layer
    timeseries = timeseries_dummy(ind_start(0):ind_end(dimsizes(ind_end)-1),0)
    dfs = dfs_dummy(ind_start(0):ind_end(dimsizes(ind_end)-1))

    ; To calculate stats of n_obs 
    ;n_values = timeseries_dummy(ind_start(0):ind_end(dimsizes(ind_end)-1),8)
    ;nstats = stat_dispersion(n_values, True)


  ; -------------------------------
  ; detrend
  ; -------------------------------
    ystart = 2002
    ymid = (tofloat(y_end-2002))/2
    detrend_val = timeseries - (global_trend*(mopitt_yfrac-ystart-ymid-0.5))
    timeseries = (/detrend_val/)

  ;--------------------------------------------
  ; Seasonal analysis
  ;--------------------------------------------
  mopitt_seas_avg_a = calc_seasonal_cycle(timeseries,time,2002,2011, 12)
  mopitt_seas_avg_b = calc_seasonal_cycle(timeseries,time,2012,y_end, 12)

  print("2002-2011 AVG --- STD DEV")
  print(mopitt_seas_avg_a(:,0)/scale + " --- " + mopitt_seas_avg_a(:,1)/scale)
  print("2012-"+y_end+" AVG --- STD DEV")
  print(mopitt_seas_avg_b(:,0)/scale + " --- " + mopitt_seas_avg_b(:,1)/scale)
  ;exit

  ;--------------------------------------------
  ; t-test
  ;--------------------------------------------
  var_a = mopitt_seas_avg_a(:,1)/scale*mopitt_seas_avg_a(:,1)/scale
  var_b = mopitt_seas_avg_b(:,1)/scale*mopitt_seas_avg_b(:,1)/scale
  num_a = mopitt_seas_avg_a(:,0)*0
  num_a = 10
  ; Note different num_a for August
  num_a(7)=9
  num_b = mopitt_seas_avg_b(:,0)*0
  num_b = 7
  df = num_a+num_b-2
  pooled_sd = sqrt((var_a+var_b)/2)
  pooled_sd_2 = sqrt(((num_a-1)*var_a+(num_b-1)*var_b)/(num_a+num_b-2))
  mean1 = mopitt_seas_avg_a(:,0)/scale
  mean2 = mopitt_seas_avg_b(:,0)/scale
  cohen_d = abs(mean1-mean2)/pooled_sd
  t_test_result = ttest(mean1,var_a,num_a,mean2,var_b,num_b,False,True)
  print(t_test_result)

  ;---- 95 % t= 2.447, 90% t = 1.934
  CI_lb = (mean1-mean2)-1.934*sqrt(var_a/num_a+var_b/num_b)
  CI_ub = (mean1-mean2)+1.934*sqrt(var_a/num_a+var_b/num_b)
  mean_diff = (mopitt_seas_avg_a(:,0)-mopitt_seas_avg_b(:,0))/scale

  ;print(mopitt_seas_avg_a(:,1)/scale)
  ;print(mopitt_seas_avg_b(:,1)/scale)

  print("*****************************")
  print("August t-test results:")
  print("p-value =   " + sprintf("%9.3f",t_test_result(0,7)) )
  print("t-value =   " + sprintf("%9.3f", t_test_result(1,7)) )
  print("df =        " + df(7))
  print("Cohen's d = " + sprintf("%9.3f", cohen_d(7)) )
  print("95%CI =     " + sprintf("%9.3f", CI_lb(7))+ ";" +sprintf("%9.3f", CI_ub(7)) )
  print("Mean diff = " + sprintf("%9.3f", mean_diff(7)) )
  print("*****************************")
  ;print("All t-test results:")
  ;print("95%CI =     " + sprintf("%9.3f",t_test_result(0,:)) )

  ;exit

sig_a = where(t_test_result(0,:).lt.(0.05), mopitt_seas_avg_a(:,0), mopitt_seas_avg_a@_FillValue)
sig_b = where(t_test_result(0,:).lt.(0.05), mopitt_seas_avg_b(:,0), mopitt_seas_avg_b@_FillValue)
sig_a01 = where(t_test_result(0,:).lt.(0.01), mopitt_seas_avg_a(:,0), mopitt_seas_avg_a@_FillValue)
sig_b01 = where(t_test_result(0,:).lt.(0.01), mopitt_seas_avg_b(:,0), mopitt_seas_avg_b@_FillValue)

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
      res@tmYLTickStartF  = miny
      res@tmYLTickEndF    = maxy
      res@tmYLTickSpacingF= 20
      ;res@tmYLTickStartF  = 160
      ;res@tmYLTickEndF    = 220

    res@trXMinF             = -0.5
    res@trXMaxF             = 11.5

    res@tmXBMode              = "Explicit"
      res@tmXBValues           = mopitt_seas_avg_a&month(1::2)
      res@tmXBLabels           = (/"Feb", "Apr", "Jun",\
                                   "Aug", "Oct", "Dec"/)
      ;res@tmXBLabels           = (/"Jan", "Feb", "Mar", "Apr", "May", "Jun",\
      ;                             "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"/)
      res@tiXAxisString          = "Year"
      ;res@tmXBLabelsOn = False
      res@tmXBLabelFontHeightF = 0.048
      res@tmYLLabelFontHeightF = 0.048
      res@tiYAxisFontHeightF  = 0.048
      res@tiXAxisFontHeightF  = 0.048

    res@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
    res@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticks inward
    res@tmXBMajorOutwardLengthF = 0.0               
    res@tmXBMinorOutwardLengthF = 0.0  

    res@tmYROn = False            ; Turn off right tickmarks.
    res@tmXTOn = False            ; Turn off top tickmarks.    

    res@vpWidthF             = 1.0
    res@vpHeightF            = 0.5          ; set height of plot
    ;res@vpYF                 = 1.0

    res@tiMainString        = ""
    res@tiXAxisString       = "Month"
    ;res@tiYAxisString       = "X~B2~CO~N~    (ppb)"
    res@tiYAxisString       = y_title
    ;res@tiYAxisString       = "DFS"
    ;res@gsnLeftString       = location + " surface"
    res@gsnLeftString       = plot_title                 ; change left 

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
    res@xyDashPatterns      =0
  plot1 = gsn_csm_xy(wks,mopitt_seas_avg_a&month,mopitt_seas_avg_a(:,0)/scale,res)


    res@xyExplicitLegendLabels   = "2012-"+y_end
    res@pmLegendOrthogonalPosF   = -1.53
    res@lgLabelFontColor       = "blue3"
    res@gsnLeftString       = ""
    res@xyLineColor         = "blue3"
    res@xyDashPatterns      =1
  plot2 = gsn_csm_xy(wks,mopitt_seas_avg_b&month,mopitt_seas_avg_b(:,0)/scale,res)
  overlay(plot1, plot2)

    res@pmLegendDisplayMode      = "Never"        ; turn on legend

  ; -------------------------------
  ; add std dev
  ; -------------------------------

  ; -------------------------------
  ; Define a polygon for std-dev
  ; -------------------------------
   month_time = mopitt_seas_avg_a&month
   ya1    = new( (/2*dimsizes(month_time)/), float )
   yb1    = new( (/2*dimsizes(month_time)/), float )
   do k=0,dimsizes(month_time)-1
     ya1(k) = month_time(k)
     yb1(k) = tofloat((mopitt_seas_avg_a(k,0)+mopitt_seas_avg_a(k,1))/scale)
     yb1(2*dimsizes(month_time)-1-k) = tofloat((mopitt_seas_avg_a(k,0)-mopitt_seas_avg_a(k,1))/scale)
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
     yb2(k) = tofloat((mopitt_seas_avg_b(k,0)+mopitt_seas_avg_b(k,1))/scale)
     yb2(2*dimsizes(month_time)-1-k) = tofloat((mopitt_seas_avg_b(k,0)-mopitt_seas_avg_b(k,1))/scale)
     ya2(2*dimsizes(month_time)-1-k) = month_time(k)
   end do
    gsres                          = True                ; poly res
    gsres@gsFillColor            = "skyblue"                ; color chosen
    gsres@gsFillOpacityF         = 0.3                 
    dum2a = gsn_add_polygon(wks,plot2,ya2,yb2,gsres)

  ; -------------------------------
  ; plot significant points
  ; -------------------------------
    res@xyMarkLineModes = "Markers"
    res@xyMarkerColors = "black"
    res@xyMarkers = 16
    res@xyMarkerSizeF    = 25.5
  plot1s = gsn_csm_xy(wks,mopitt_seas_avg_a&month,sig_a/scale,res)
    res@xyMarkerSizeF    = 25.5
  ;plot1s01 = gsn_csm_xy(wks,mopitt_seas_avg_a&month,sig_a01/scale,res)
  overlay(plot1, plot1s)
  ;overlay(plot1, plot1s01)

    res@xyMarkerColors = "blue3"
    res@xyMarkerSizeF    = 25.5
  plot2s = gsn_csm_xy(wks,mopitt_seas_avg_b&month,sig_b/scale,res)
    res@xyMarkerSizeF    = 25.5
  ;plot2s01 = gsn_csm_xy(wks,mopitt_seas_avg_a&month,sig_b01/scale,res)
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
  ;  pm01res = True
  ;  pm01res@gsMarkerIndex = 16
  ;  pm01res@gsMarkerSizeF = 0.020
  ;  pm01res@gsMarkerColor := "black"
   ; gsn_polymarker_ndc (wks,0.15,0.72,pm01res)
   ; gsn_text_ndc(wks,"p < 0.01",0.214,0.72,txres)


  ; -------------------------------
  ; ADD CUSTOM TEXT
  ; -------------------------------
  ;drawNDCGrid(wks)
  ; STATIONS
  txres                   = True                ; polyline mods desired
  txres@txFontHeightF = 0.028                   ; font smaller. default big
  txres@txJust = "CenterLeft"                   ; font smaller. default big
  txres@txFontColor = "black"
  ;gsn_text_ndc(wks,"MOPITT 2002-2011",0.60,0.69,txres)
  txres@txFontColor = "blue3"
  ;gsn_text_ndc(wks,"MOPITT 2012-"+y_end,0.60,0.66,txres)

  ; -------------------------------
  ; DRAW IT OUT
  ; -------------------------------
  draw(plot1)
  frame(wks)

end if



end

