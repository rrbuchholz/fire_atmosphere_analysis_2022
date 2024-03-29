;=======================================================;
;  regional-subset-mopitt-stats.ncl
;=======================================================;
;
; Concepts illustrated:
;   - For use with MOPITT CO version 6,7 and 8 level 2 products
;   - Loading CO data from MOPITT .he5 file
;          - select daytime retrievals
;          - select for spatial region
;          - select land or ocean only pixels
;          - determine, mean, median, 5th, 
;            25th, 75th, 95th percentiles
;            in month average values
;   - Write data out to netCDF
;
;   Files needed:
;         -MOPITT level 2, daily hdf files
;
;   To run type:
;         ncl regional-subset-mopitt.ncl 
;                                  RRB January 8, 2015
;=======================================================;

begin

;--------------------------------------------
; user input
;--------------------------------------------
  ystart = 2002
  yend = 2018
  vnumber = "V8"
  retrieval = "J"
  version = vnumber+retrieval
  infolder = "/MOPITT/"+version+"/Archive/L2/"

  ;------------
  ; select location
  ;------------
  location = "PNW"    ; region

  ;------------
  ; select daytime 
  ; retrievals
  ;------------
  ; For SZA < 80
  dayflag = 80

  ;------------
  ; select surface type 
  ;------------
  surface = "land"

  ;------------
  ; toggles
  ;------------
  netCDF         = True
  COLAVG_VMR     = True

;--------------------------------------------
; end user input
;--------------------------------------------
;
;--------------------------------------------
; set up
;--------------------------------------------
  ; -------------------------------
  ; MOPITT
  ; -------------------------------
  ; names of data structures
  ; determined from an ncl_filedump
  suff           = "_MOP02"
  tracer         = "RetrievedCOTotalColumn"+suff
  longitude      = "Longitude"+suff
  latitude       = "Latitude"+suff
  date           = "Time"+suff
  deg_freedom    = "DegreesofFreedomforSignal"+suff
  tcol_apriori   = "APrioriCOTotalColumn"+suff
  vmr_aprior     = "APrioriCOMixingRatioProfile"+suff
  vmr_surface    = "APrioriCOSurfaceMixingRatio"+suff
  vmr_aposterior = "RetrievedCOMixingRatioProfile"+suff
  apost_surf     = "RetrievedCOSurfaceMixingRatio"+suff
  avker          = "TotalColumnAveragingKernel"+suff
  vmr_avker      = "RetrievalAveragingKernelMatrix"+suff
  pressure_surf  = "SurfacePressure"+suff
  diag           = "RetrievedCOTotalColumnDiagnostics"+suff
  dryair         = "DryAirColumn"+suff

  ; used for selecting, but not written out:
  solarza        = "SolarZenithAngle"+suff
  surface_index  = "SurfaceIndex"+suff
  diag_2         = "RetrievalAnomalyDiagnostic"+suff
  swath_index    = "SwathIndex"+suff
  snr_calc       = "Level1RadiancesandErrors"+suff


   region_names = (/"PNW", "CentralUSA", "EastUSA", "Colorado", "Northeast"/)

   ;                   minlat, maxlat, minlon, maxlon
   region_select = (/(/ 38.,  57.,     -127.,     -110./),\ ;PNW
                     (/ 35.,  49.,     -110.,     -95./),\  ;Central Mountain USA
                     (/ 32.5,  42.,     -95.,     -75./),\  ;East USA
                     (/37.,    41.,    -109.,    -102./),\  ;Colorado
                     (/ 41.,  53.,     -95.,     -74./) /)  ;Northeast


   ind_region = ind(region_names.eq.location)
   topboundary     = region_select(ind_region,1)
   bottomboundary  = region_select(ind_region,0)
   rightboundary   = region_select(ind_region,3)
   leftboundary    = region_select(ind_region,2)
  
  ;Carriage return
  cr = tochar(10)

  ; -------------------------------
  ; SELECTION SETUP
  ; -------------------------------
  if (surface.eq."land") then
     s_id = 1
  else if (surface.eq."water") then
     s_id = 0
  end if
  end if

do y = ystart,yend ; loop through years
  year = y
  ; -------------------------------
  ; SETUP OUTPUT
  ; -------------------------------
  if (netCDF) then
     outfolder = "~/MOPITT_subset/regions_v8/"
     outname      = outfolder+location+"_"+version+"subset_"+year+"monthavg.nc"
     if (COLAVG_VMR) then
        outname      = outfolder+location+"_"+version+"subset_"+year+"monthavg_VMR.nc"
     end if
  end if

  ; -------------------------------
  ; SETUP ARRAYS 
  ; -------------------------------
  yearmonth = yyyymm_time(year,year, "integer")
  yearmonth_string = tostring(yearmonth)
  month_stats = new((/dimsizes(yearmonth_string),12/), float)
  month_5_95 = new((/dimsizes(yearmonth_string),2/), float)
  month_tcol = new((/dimsizes(yearmonth_string)/), float)
  month_err = new((/dimsizes(yearmonth_string)/), float)
  month_rande = new((/dimsizes(yearmonth_string)/), float)
  month_smoothe = new((/dimsizes(yearmonth_string)/), float)
  month_dofs = new((/dimsizes(yearmonth_string)/), float)
  month_apriori = new((/dimsizes(yearmonth_string)/), float)
  month_psurf = new((/dimsizes(yearmonth_string)/), float)
  month_dryair = new((/dimsizes(yearmonth_string)/), float)
  month_ret_surf = new((/dimsizes(yearmonth_string),2/), float)

  month_prior = new((/dimsizes(yearmonth_string), 10/), float)
  month_prof = new((/dimsizes(yearmonth_string), 10/), float)
  month_colak = new((/dimsizes(yearmonth_string), 10/), float)

  month_vmrak = new((/dimsizes(yearmonth_string), 10, 10/), float)

;--------------------------------------------
; load file and extract
;--------------------------------------------
  ;------------
  ; loop over
  ; months
  ;------------
  do m=0,11
    ; counter for counting number of files with selected values
    dcount = 0
    print("==================================")
    print("Collecting: "+ yearmonth_string(m))
    ; -------------------------------
    ; MOPITT
    ; -------------------------------
    all_files := systemfunc ("ls "+infolder+yearmonth_string(m)+"/*/MOP02"+retrieval+"-"+year+"*.he5")
    ; number of files to be read
    ndates = dimsizes(all_files)

    if (all(ismissing(all_files))) then
      print("No MOPITT files for date: "+cd_calendar(time(i),2))
      continue  ; skip all the next bit and go to next file
    end if

    ;------------
    ; loop over
    ; files in month
    ;------------
    do i=0,ndates-1
      ; read data
      fin = addfile(all_files(i), "r")            ; read in
      sza          := fin->$solarza$
      surf_type    := fin->$surface_index$
      pixel_num    := fin->$swath_index$(:,0)
      sn_5A        := fin->$snr_calc$(:,3,:)
      sn_6A        := fin->$snr_calc$(:,9,:)
        snr_5A := sn_5A(:,0)/sn_5A(:,1)
        snr_6A := sn_6A(:,0)/sn_6A(:,1)
      lat          := fin->$latitude$
      lon          := fin->$longitude$
      diag_anom    := fin->$diag_2$
      time         := fin->$date$
      dofs         := fin->$deg_freedom$
        ;units value presumes use of TAI93 (International Atomic Time) format
        time@units = "seconds since 1993-1-1 00:00:00"


      ;--------------------------------------------
      ; Select region to bin (daytime retrievals)
      ;--------------------------------------------
      region := ind(sza.le.dayflag.and.\
                surf_type.eq.s_id.and.\
                pixel_num.ne.3.and.\
                snr_5A.ge.1000.and.\
                ;snr_6A.ge.400.and.\
                lat.le.topboundary.and.\
                lat.ge.bottomboundary.and.\
                lon.ge.leftboundary.and.\
                lon.le.rightboundary.and.\
                diag_anom(:,0).ne.1.and.\
                diag_anom(:,1).ne.1.and.\
                diag_anom(:,2).ne.1.and.\
                diag_anom(:,3).ne.1.and.\
                diag_anom(:,4).ne.1)

      if (all(ismissing(region))) then
        print("No MOPITT data at " + location +\
              " for file " + i +", date: "+cd_calendar(time(i),2))
        continue  ; skip all the next bit and go to next file
      end if

      if (dimsizes(region).eq.1) then
        print("only one value from MOPITT at " + location +\
              " for file " + i +", date: "+cd_calendar(time(i),2))
        continue  ; skip all the next bit and go to next file
      end if

      tempcount = dcount
      dcount := tempcount+1
      print("***Collecting MOPITT data at " + location +\
              " for file " + i +", date: "+cd_calendar(time(i),2))

      ;--------------------------------------------
      ; Load remaining variables
      ;--------------------------------------------
      tgas          := fin->$tracer$(region,:)
        region_tcol := tgas(:,0)
        region_err  := tgas(:,1)
        region_two_diag  := fin->$diag$(region,:)
          region_smooth_e:= max(abs(region_two_diag(:,0)))
          region_rand_e  := max(abs(region_two_diag(:,1)))
        region_dofs       := fin->$deg_freedom$(region)
        region_apriori    := fin->$tcol_apriori$(region)
        region_psurf      := fin->$pressure_surf$(region)
        region_dryair     := fin->$dryair$(region)

     ; Convert to column average VMR
     if (COLAVG_VMR) then
       col_vmr      := region_tcol/region_dryair
       col_vmr_err  := region_err/region_dryair
       ap_vmr       := region_apriori/region_dryair
       region_tcol  := col_vmr
       region_err   := col_vmr_err
       region_apriori := ap_vmr
     end if 

      ; a priori and retrieved profile have surface layer separate to rest
      vmr_dummy    := fin->$vmr_aprior$(region,:,0)
      vmr_surf     := fin->$vmr_surface$(region,0)
      region_prior := new((/dimsizes(region_tcol),10/),float,-9999)
         region_prior(:,1:9)  = vmr_dummy(:,:)
      ret_dummy    := fin->$vmr_aposterior$(region,:,0)
      ret_surf     := fin->$apost_surf$(region,0)
      region_prof  := new((/dimsizes(region_tcol),10/),float,-9999)
         region_prof(:,1:9)  = ret_dummy(:,:)
      region_ak    := fin->$avker$(region,:)
      region_vmr_ak:= fin->$vmr_avker$(region,:,:)

     ;--------------------------------------------------------------
     ; Correct for where MOPITT surface pressure <900 hPa
     ;--------------------------------------------------------------
     dp := new((/dimsizes(region_psurf), 10/), float) ; Pressure difference array
     do z= 0, 9, 1
       dp(:,z) = 1000.-(z*100) - region_psurf
     end do
     ;print(region_psurf)

     do p = 0, dimsizes(region_tcol)-1
       region_prior(p,:) = where(ismissing(region_prior(p,:)),vmr_surf(p),region_prior(p,:))
       region_prof(p,:) = where(ismissing(region_prof(p,:)),ret_surf(p),region_prof(p,:))
     end do
     region_prior = where(dp.ge.100,region_prior@_FillValue,region_prior)
     region_prof = where(dp.ge.100,region_prof@_FillValue,region_prof)


      if (any(.not.ismissing(region_tcol))) then
      ;--------------------------------------------
      ; Collect MOPITT data into timeseries
      ;--------------------------------------------
          if (.not.(isdefined("mopitt_selected"))) then

            ;--- total column
            mopitt_selected  = region_tcol
            print("New data points: " +dimsizes(region_tcol)+\
                  ", Total data points "+dimsizes(mopitt_selected))

            ;--- 1D vars
            dofs_selected    = region_dofs
            apriori_selected = region_apriori
            psurf_selected   = region_psurf
            dryair_selected  = region_dryair
            err_selected     = region_err
            rande_selected   = region_rand_e
            smoothe_selected = region_smooth_e
            ret_surf_selected = ret_surf

            ;--- 2D vars
            m_apriori_selected = region_prior
            m_profile_selected = region_prof
            ak_selected        = region_ak

            ;--- 3D vars
            vmr_ak_selected    = region_vmr_ak

          else
            ;--- total column
            temp1 := mopitt_selected
            mopitt_selected := array_append_record(temp1,region_tcol,0)

            ;--- 1D vars
            temp2 := dofs_selected
            temp3 := apriori_selected
            temp4 := psurf_selected
            temp5 := dryair_selected
            temp6 := err_selected
            temp7 := rande_selected
            temp8 := smoothe_selected
            temp8a := ret_surf_selected
            dofs_selected    := array_append_record(temp2,region_dofs,0)
            apriori_selected := array_append_record(temp3,region_apriori,0)
            psurf_selected   := array_append_record(temp4,region_psurf,0)
            dryair_selected  := array_append_record(temp5,region_dryair,0)
            err_selected     := array_append_record(temp6,region_err,0)
            rande_selected   := array_append_record(temp7,region_rand_e,0)
            smoothe_selected := array_append_record(temp8,region_smooth_e,0)
            ret_surf_selected := array_append_record(temp8a,ret_surf,0)

            ;--- 2D vars
            temp9 := m_apriori_selected
            temp10 := m_profile_selected
            temp11 := ak_selected
            m_apriori_selected := array_append_record(temp9,region_prior,0)
            m_profile_selected := array_append_record(temp10,region_prof,0)
            ak_selected        := array_append_record(temp11,region_ak,0)

            ;--- 3D vars
            temp12 := vmr_ak_selected
            vmr_ak_selected := array_append_record(temp12,region_vmr_ak,0)

            print("New data points: " +dimsizes(region_tcol)+\
                  ", Total data points "+dimsizes(mopitt_selected)+\
                  ", check1 "+dimsizes(dofs_selected)+\
                  ", check2 "+dimsizes(apriori_selected))
          end if
      end if

      end do ; file loop

      ; still no data
      if (.not.(isdefined("mopitt_selected"))) then
        print("No MOPITT data at " + location +\
              " for " +", date: "+yearmonth_string(m))
        continue  ; skip all the next bit and go to next file
      end if

    ;--------------------------------------------
    ; Perform percentile calculations
    ;--------------------------------------------
      opt = False
      stats_array = stat_dispersion(mopitt_selected, opt )
      sorted_array := dim_pqsort(mopitt_selected,2) 
      ind_5th = toint(ceil(0.05*dimsizes(sorted_array)))
      ind_95th = toint(ceil(0.95*dimsizes(sorted_array)))
      if (ind_95th.ge.dimsizes(mopitt_selected))then
        ind_95th = dimsizes(mopitt_selected)-1
      end if

      ; Stats Output
      ;(0)	   ===> Robust Dispersion Statistics: RetrievedCOTotalColumn <===
      ;(0)	 [0]            Mean
      ;(0)	 [1]          StdDev
      ;(0)	 [2]             Min
      ;(0)	 [6]     LowQuartile
      ;(0)	 [8]          Median
      ;(0)	 [10]   HighQuartile
      ;(0)	 [14]            Max
      ;(0)	 [18]      #   Total
      ;(0)	 [19]      #    Used=
      ;(0)	 [20]      # Missing
      ;(0)	 [24]     Lower 5.0%=-9999   <------- n > 1000
      ;(0)	 [25]     Upper 5.0%=-9999   <------- n > 1000
    ;------------
    ; collect month
    ; stats
    ;------------
    ; tcol
    month_stats(m,:) = stats_array((/0,1,2,6,8,10,14,18,19,20,24,25/))
    month_5_95(m,0) = mopitt_selected(ind_5th)
    month_5_95(m,1) = mopitt_selected(ind_95th)
    month_tcol(m) = dim_avg_n_Wrap(mopitt_selected,0)

    ; 1D vars
    month_dofs(m) = dim_avg_n_Wrap(dofs_selected,0)
    month_apriori(m) = dim_avg_n_Wrap(apriori_selected,0)
    month_psurf(m) = dim_avg_n_Wrap(psurf_selected,0)
    month_dryair(m) = dim_avg_n_Wrap(dryair_selected,0)
    ; errors are summed in quadrature
    month_err(m) = sqrt(dim_avg_n_Wrap(err_selected^2,0))
    month_rande(m) = sqrt(dim_avg_n_Wrap(rande_selected^2,0))
    month_ret_surf(m,0) = dim_avg_n_Wrap(ret_surf_selected,0)
    month_ret_surf(m,1) = dim_stddev_n(ret_surf_selected,0)
    
    ; 2D vars
    month_prior(m,:) = dim_avg_n_Wrap(m_apriori_selected,0)
    month_prof(m,:) = dim_avg_n_Wrap(m_profile_selected,0)
    month_colak(m,:) = dim_avg_n_Wrap(ak_selected,0)

    ; 3D vars
    month_vmrak(m,:,:) = dim_avg_n_Wrap(vmr_ak_selected,0)

    records = num(mopitt_selected)
    print("Number of files checked = " + ndates)
    print("Number of files with " +location+" info for "+ yearmonth_string(m)+" = " + dcount)
    print("Number of retrievals in month values = " + records)
    
    ;------------
    ; release variable
    ;------------
      delete([/mopitt_selected,dofs_selected,apriori_selected,psurf_selected,\
              dryair_selected,err_selected,rande_selected,smoothe_selected,\
              m_apriori_selected,m_profile_selected,ak_selected,vmr_ak_selected,\
              ret_surf_selected/])


;--------------------------------------------
; Write out to a netCDF file
;--------------------------------------------
if (netCDF) then

  ; saves MOPITT subset
    print("Creating file...")
    fout = addfile(outname, "c")
    ;------------
    ; File 
    ; Attributes
    ;------------
    setfileoption(fout, "DefineMode",True)
    fAtt                      = True
      fAtt@title              = "Subset of month average MOPITT CO for " +location+"_"+year
      fAtt@source             = "MOPITT level 2, "+version+" retrievals"
      fAtt@creation_date      = systemfunc ("date")
      fAtt@conventrions       = "CF"

     fileattdef(fout, fAtt)
     setfileoption(fout, "DefineMode",False)


  ;------------
  ; Variables
  ;------------
   month_tcol!0              = "time"
   month_tcol&time           = yearmonth
   month_tcol@average_op_ncl = "monthly spatial average over "+location+": "+\
                                    topboundary+" to "+bottomboundary+ " Lat, "+\
                                    leftboundary+" to "+rightboundary+" Lon"
   month_stats!0               = "time"
   month_stats&time            = yearmonth
   month_stats@average_op_ncl  = "monthly spatial statistics over "+location+": "+\
                                    topboundary+" to "+bottomboundary+ " Lat, "+\
                                    leftboundary+" to "+rightboundary+" Lon"
   month_stats@definitions ="[0] Mean, [1] StdDev, [2] Min, [3] LowQuartile, [4] Median, "+cr+\
       "[5] HighQuartile, [6] Max, [7] Total n, [8] Used n, [9] Missing, "+cr+\
       "[10] Lower 5.0% (for n > 1000), [11] Upper 5.0% (for n > 1000)"

   month_5_95!0               = "time"
   month_5_95&time            = yearmonth
   month_5_95@average_op_ncl  = "monthly 5th and 95th percentiles over "+location+": "+\
                                    topboundary+" to "+bottomboundary+ " Lat, "+\
                                    leftboundary+" to "+rightboundary+" Lon"

   month_apriori!0               = "time"
   month_apriori&time            = yearmonth

   if (COLAVG_VMR) then
     ;--- total column
       fout->RetrievedX_CO       = month_tcol
     ;--- stats
       fout->RetrievedX_CORegionStats    = month_stats
       fout->RetrievedX_CORegion5th95th  = month_5_95
       fout->AvgAPrioriX_CO              = month_apriori
   else
     ;--- total column
       fout->RetrievedCOTotalColumn        = month_tcol
     ;--- stats
       fout->RetrievedCOtcolRegionStats    = month_stats
       fout->RetrievedCOtcolRegion5th95th  = month_5_95
       fout->AvgAPrioriCOTotalColumn      = month_apriori
   end if
   ;--- 1D vars
   month_dofs!0               = "time"
   month_dofs&time            = yearmonth
     fout->AvgDegreesofFreedomforSignal = month_dofs
   month_psurf!0               = "time"
   month_psurf&time            = yearmonth
     fout->AvgSurfacePressure           = month_psurf
   month_dryair!0               = "time"
   month_dryair&time            = yearmonth
     fout->AvgDryAirColumn              = month_dryair
   month_err!0               = "time"
   month_err&time            = yearmonth
     fout->AvgError                     = month_err
   month_rande!0               = "time"
   month_rande&time            = yearmonth
     fout->AvgRandomError               = month_rande
   month_smoothe!0               = "time"
   month_smoothe&time            = yearmonth
     fout->AvgSmoothingError            = month_smoothe
   month_ret_surf!0               = "time"
   month_ret_surf&time            = yearmonth
   month_ret_surf@definitions     = "[:,0] Mean [:,1] StdDev"
     fout->AvgSurfaceVMR          = month_ret_surf

   ;--- 2D vars
   month_prior!0               = "time"
   month_prior&time            = yearmonth
   month_prior!1               = "lev"
     fout->AvgAPrioriCOMixingRatioProfile = month_prior
   month_prof!0               = "time"
   month_prof&time            = yearmonth
   month_prof!1               = "lev"
     fout->AvgRetrievedCOMixingRatioProfile = month_prof
   month_colak!0               = "time"
   month_colak&time            = yearmonth
   month_colak!1               = "lev"
     fout->AvgTotalColumnAveragingKernel = month_colak

   ;--- 3D vars
   month_vmrak!0               = "time"
   month_vmrak&time            = yearmonth
     fout->AvgRetrievalAveragingKernelMatrix = month_vmrak


end if ;NETCDF

end do ; year loop

end
