; This is an NCL script file which plots residuals
; between the EwE simulated and EwE-F simulated
; Generic 37 model absolute biomass values
; in a box-whisker plot fashion.
;************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;************************************
begin

 ewefort = readAsciiTable("EwE-F_Generic_37_Biomass.csv", 38, "float", 1)
 ewe6    = readAsciiTable("EwE-6.6_Generic_37_Biomass.csv", 38, "float", 10)

 residuals = ewe6(:, 1:37) - ewefort(1:120, 1:37)
 residualsT = transpose(residuals)

;************************************
;create plot
;************************************

 wks_type = "pdf"
 wks_type@wkPaperSize = "A4"
 wks_type@wkOrientation = "landscape"
 wks = gsn_open_wks(wks_type, "Generic_37_residuals") ; open a ps  file

 xdata = ispan(1, 37, 1)
 iarr  = new((/37, 5/), float)            ;  fill with minimum, 25th percentile, median, 75th percentile, maximum of each timeseries

 do i = 0, 36
     dimt = dimsizes(residualsT(i, :))
     x25  = round(.25 * dimt, 3) - 1     ; -1 to account for NCL indexing starting
     x75  = round(.75 * dimt, 3) - 1     ; at 0            
     qsort(residualsT(i, :))             ; sort the data
     iarr(i, :) = (/min(residualsT(i, :)), residualsT(i, x25), dim_median(residualsT(i, :)), residualsT(i, x75), max(residualsT(i, :))/)  
 end do
 
 res = True

 res@gsnMaximize = True        ; Maximize plot in frame; aspect ratio will be preserved.
 res@vpWidthF  = 0.8           ; Set width and height of plot.
 res@vpHeightF = 0.3

 res@tiMainString  = "Generic 37 Model"
 res@tiMainFontHeightF = 0.01
 res@tiYAxisString = "~F8~d"
 res@tmXBLabels    = (/"Baleen whales", "Toothed whales", "Seals", "Birds", "Sharks, large", "Sharks, small medium", "Rays, large", "Rays, small medium", \
                       "Pelagics, large", "Pelagics, medium", "Pelagics, small, carniv.", "Pelagics, small, herbiv.", "Benthopelagics, large", "Benthopelagics, small medium", \
                       "Demersals, large", "Demersals, medium", "Demersals, small", "Reeffish, large", "Reeffish, medium", "Flatfish, large", "Flatfish, small medium", \
                       "Reeffish, small", "Bathypelagics", "Bathydemersals", "Jellyfish", "Cephalopods", "Shrimps", "Lobsters, crabs", "Macrobenthos", "Meiobenthos", \
                       "Corals", "Softcorals, sponges, etc", "Krill", "Zooplankton, other", "Phytoplankton", "Benthic plants","Detritus"/)

 res@tmXBLabelFontHeightF    = 0.009        ; Increase font height
 res@tmXBLabelAngleF         = -90
 res@tmYLLabelFontHeightF    = 0.013        ; Increase font height
 res@tmYLFormat              = "!;*^se"     ; Exponential notation for Y axis
 res@trYMaxF                 = 0.0001
 res@trYMinF                 = -0.0001
 res@trXMaxF                 = 38
 res@trXMinF                 = 0
 res@tmYLMajorOutwardLengthF = 0.02

;**********************************************
; resources for polylines that draws the boxes
;**********************************************  
;  llres                   = True			
;  llres@gsLineThicknessF  = 0.5                 ; line thickness 
;**********************************************
; resources that control color and width of boxes
;**********************************************  
;  opti          = True			
;  opti@boxWidth = 0.25				 ; Width of box (x units)
;***********************************************

 plot = boxplot(wks, xdata, iarr, False, res, False)	
 draw(wks)                                    ; boxplot does not call these
 frame(wks)                                   ; for you

end

