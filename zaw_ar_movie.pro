;----------------------------------------------------------------------------------------------------------
; NAME:           
;             zaw_ar_movie
;
; PURPOSE:        
;             Inputs a filename and plots the magnetograms with AR overlays
;
; CALLING SEQUENCE:
;             zaw_ar_movie, filename    
;
; INPUTS:
;             filename                          - file name of the magnetograms to plot from
;                 
;
; OPTIONAL INPUTS:
;             start_date                        - where the plotting begins
;             end_date                          - where the plotting ends
;       
;             /prt                              - does something
;                  
; OUTPUT:
;             [none]
;
;-------------------------------------------------------------------------------------
PRO zaw_ar_movie, filename=filename, instr, start_date = start_date, end_date = end_date, prt = prt

;Set initial parameters
IF( keyword_set( filename ) ) THEN file = filename $
ELSE file = 'Ar_id_Save0.sav'

restore, file
ARs = ARs[1:n_elements(ARs)-1]

if keyword_set(start_date) then begin
	min_date = long( mdi_datestr( start_date ) )
endif else begin
	min_date = min(ARs.mdi_i)
endelse

if keyword_set(end_date) then begin
	max_date = long( mdi_datestr( end_date ) )
endif else begin
max_date = max(ARs.mdi_i)
endelse

;Window Initialization
pls = 500; Magnetogram size
pad = 40;  Padding
plxy = 0; Bottom space
plxx = 0; Right space

;First reference day for keeping track time easily
if instr eq 1 then DayOff = julday(1,1,1970); KPVT 512
if instr eq 2 then DayOff = julday(1,1,1990); KPVT SPMG
if instr eq 3 then DayOff = julday(1,1,1993); MDI
if instr eq 4 then DayOff = julday(1,1,2009); HMI

;Main Loop
for i = min_date, max_date do begin
   
  lw_msw = 1
  n_itr = 0
  mdi_il = i
  if lw_msw ne 0 then begin
    print, 'Looking for valid Reference Window magnetogram' 
    repeat begin
        ;Reading files 
        caldat, mdi_il+DayOff, Month, Day, Year
        datel = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)
        print, datel
        mgl = amj_file_read( datel, hdrl, instr )
    
        sl = size(mgl)
    
      mdi_il = mdi_il + lw_msw
      n_itr = n_itr+1
      
    endrep until ( (sl[2] eq 8) or (n_itr gt 365) )
  
    mdi_il = mdi_il - lw_msw
    lw_msw = 0

  endif

  i = mdi_il
  
  sz=size(mgl.img)
  d_zoom = pls/float(sz[1])
  
  ;window creations
  set_plot,'X'
  window,0,xsize=pls+4.0*pad+plxx,ysize=pls+2.0*pad+plxy
  
  ;Exploration Window
  ;AR overlay
  if keyword_set(prt) then begin
	amj_mgplot, mgl.img, mdi_il, ARs = ARs, d_xsize = pls, d_ysize = pls, max = 300, pos = [2.0*pad, pad+plxy, pls+2.0*pad, pls+pad+plxy], title = datel, xrange = [min(mgl.x), max(mgl.x)] , yrange = [min(mgl.y), max(mgl.y)], sqs_nm = i-min_date+1, /prt, /shw_lbl   
  endif else begin
	amj_mgplot, mgl.img, mdi_il, ARs = ARs, d_xsize = pls, d_ysize = pls, max = 300, pos = [2.0*pad, pad+plxy, pls+2.0*pad, pls+pad+plxy], title = datel, xrange = [min(mgl.x), max(mgl.x)] , yrange = [min(mgl.y), max(mgl.y)], /shw_lbl
  endelse
  
  Wait, 1.0 ; pause between each frame in the movie
  
endfor
  
return
END