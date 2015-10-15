;----------------------------------------------------------------------------------------------------------

PRO mdd_ar_kpvt_movie, filename = filename, start_date = start_date, end_date = end_date, prt = prt


;initializing run and variables
;-------------------------------------------------------------------------------------

IF( keyword_set( filename ) ) THEN file = filename $
ELSE file = 'ARs0.sav'
 
restore, file

ARs = ARs[1:n_elements(ARs)-1]
DayOff = julday(1,1,1970);  First reference day for keeping track time easily


if keyword_set(start_date) then begin
;  minimum = long( mdi_datestr( start_date ) )
  minimum = julday(strmid(start_date,5,2),strmid(start_date,8,2),strmid(start_date,0,4))-DayOff
endif else begin
	minimum = min(ARs.mdi_i)
endelse

if keyword_set(end_date) then begin
;  maximum = long( mdi_datestr( end_date ) )
  maximum = julday(strmid(end_date,5,2),strmid(end_date,8,2),strmid(end_date,0,4))-DayOff
endif else begin
	maximum = max(ARs.mdi_i)
endelse



for i = minimum, maximum do begin
   
   mdi_il = i;

  ;Reading files
  ;MDI
;  datel = mdi_datestr( string( mdi_il ), /inverse ) 
;  mgl = mdi_archive( datel, i=1, win=[-1800, -1800, 1800, 1800], hdrl )

  ;KPVT
  print, 'Looking for valid Reference Window magnetogram' 
  repeat begin
    ;Reading files 
    caldat, mdi_il+DayOff, Month, Day, Year
    datel = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)
    print, datel
    mgl = amj_kpvt_archive( datel, hdrl )
  
    sl = size(mgl)
  
    mdi_il = mdi_il + lw_msw
    n_itr = n_itr+1
    
  endrep until ( (sl[2] eq 8) or (n_itr gt 365) )
  
  mdi_il = mdi_il - 1
  i = mdi_il


   
  ;Definition of window parameters and window initialization
  pls = 500; Magnetogram size
  pad = 40;  Padding
  plxy = 0; Bottom space
  plxx = 0; Right space
  
  sz=size(mgl.img)
  d_zoom = pls/float(sz[1])
  
  ;window creations
  set_plot,'X'
  window,0,xsize=pls+4.0*pad+plxx,ysize=pls+2.0*pad+plxy
  
  ;Exploration Window
  ;AR overlay
  if keyword_set(prt) then begin
	amj_mgplot, mgl.img, mdi_il, ARs = ARs, d_xsize = pls, d_ysize = pls, max = 300, pos = [2.0*pad, pad+plxy, pls+2.0*pad, pls+pad+plxy], title = datel, xrange = [min(mgl.x), max(mgl.x)] , yrange = [min(mgl.y), max(mgl.y)], sqs_nm = i-minimum+1, /prt, /shw_lbl   
  endif else begin
	amj_mgplot, mgl.img, mdi_il, ARs = ARs, d_xsize = pls, d_ysize = pls, max = 300, pos = [2.0*pad, pad+plxy, pls+2.0*pad, pls+pad+plxy], title = datel, xrange = [min(mgl.x), max(mgl.x)] , yrange = [min(mgl.y), max(mgl.y)], /shw_lbl
  endelse
  
  Wait, 1.0 ; pause between each frame in the movie
  
endfor
  
return
END


;----------------------------------------------------------------------------------------------------------

 
