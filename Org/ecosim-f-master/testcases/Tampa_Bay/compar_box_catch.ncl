; This is an NCL script file which plots residuals
; between the EwE simulated and EwE-F simulated
; Tampa Bay model absolute biomass values
; in a box-whisker plot fashion.
;************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;************************************
begin

 ewefort = readAsciiTable("EwE-F_Tampa_Bay_Yield.csv", 53, "float", 1)
 ewe6    = readAsciiTable("EwE-6.6_Tampa_Bay_Yield.csv", 53, "float", 10)

 residuals = ewe6(:, 1:52) - ewefort(1:612, 1:52)
 residualsT = transpose(residuals)

;************************************
;create plot
;************************************

 wks_type = "pdf"
 wks_type@wkPaperSize = "A4"
 wks_type@wkOrientation = "landscape"
 wks = gsn_open_wks(wks_type, "Tampa_Bay_catches_residuals") ; open a ps  file

 xdata = ispan(1, 52, 1)
 iarr  = new((/52, 5/), float)            ;  fill with minimum, 25th percentile, median, 75th percentile, maximum of each timeseries

 do i = 0, 51
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

 res@tiMainString  = "Tampa Bay Model Catches"
 res@tiMainFontHeightF = 0.01
 res@tiYAxisString = "~F8~d"
 res@tmXBLabels    = (/"0-12 Snook", "3-12 Snook", "12-48 Snook", "48-90 Snook", "90+ Snook", "0-3 Red Drum", \
               "3-8 Red Drum", "8-18 Red Drum", "18-36 Red Drum", "36+ Red Drum", "0-3 Sea Trout", \
               "3-18 Sea Trout", "18+ Sea Trout", "0-3 Sand Trout", "3-12 Sand Trout", "12+ Sand Trout", \
               "0-6 Mullet", "6-18 Mullet", "18+ Mullet", "Mackrel 0-3", "Mackrel 3+", "Ladyfish 0-10", \
               "Ladyfish 10+", "Jacks", "Bay Anchovy", "Pin Fish", "Spot", "Silver Perch", "Scaled Sardine", \
               "Mojarra", "Threadfin Herring", "Manhaden", "Menidia (silverside)", "Catfish", "Bumper", "Caridan Shrimp", \
               "Shrimp", "Stone Crab", "Blue Crab", "Cyprinodontids", "Poecilids", "Pigfish", "Gobies", "Rays", \
               "Benthic Invertebrates", "Macro Zooplankton", "Micro Zoolplankton", "Infauna", "Attached Microalgae", \
               "Sea Grass", "Phytoplankton", "Detritus"/)

 res@tmXBLabelFontHeightF    = 0.009        ; Increase font height
 res@tmXBLabelAngleF         = -90
 res@tmYLLabelFontHeightF    = 0.013        ; Increase font height
 res@tmYLFormat              = "!;*^se"     ; Exponential notation for Y axis
 res@trYMaxF                 = 0.0001
 res@trYMinF                 = -0.0001
 res@trXMaxF                 = 53
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

