
FUNCTION button_choice, pls, pad, plxy, plxx, dbl_sw

ftsz = 14

;Defining control (right) buttons
;labs = ['RE-DETECT', 'RE-LABEL', 'AGGREGATE', 'CREATE', 'DELETE', 'QUIT' ]
labs = ['PRINT SCREEN', 'SYNCHRONIZE','CONTRAST -','CONTRAST +', 'QUIT' ]

nl = n_elements( labs )

x0c = (1.0 + dbl_sw)*pls+(2.0*(1.0 + dbl_sw) + 0.5)*pad
x1c = (1.0 + dbl_sw)*pls+(2.0*(1.0 + dbl_sw) - 0.5)*pad+plxx
ybc = pad + plxy + float(findgen(nl+1))/float(nl)*pls


;Plotting control (right) buttons
plots, x0c*[1, 1],  [ybc[0], ybc[nl]], /device
plots, x1c*[1, 1],  [ybc[0], ybc[nl]], /device

plots, [ x1c,  x0c ], ybc[0]*[1, 1], /device
FOR i=0, nl-1 DO BEGIN
  y1 = ybc[i+1]
  y0 = 0.5*( ybc[i] + ybc[i+1] ) - ftsz/2
  xyouts, 0.5*(x1c+x0c), y0, labs[i], /device, align=0.5, charsize=2
  plots, [ x1c,  x0c ], y1*[1, 1], /device
ENDFOR



;Defining (bottom left) buttons
labs2 = ['< x10', '<', 'JUMP TO','>',  '> x10']
n2 = n_elements( labs2 )

xbr = 2.0*pad + float(findgen(n2+1))/float(n2)*pls

y0r = pad/4
y1r = plxy - pad*3/4

;Plotting revision (bottom left) buttons
plots, [xbr[0], xbr[n2]], y0r*[1, 1], /device
plots, [xbr[0], xbr[n2]], y1r*[1, 1], /device

plots, xbr[0]*[1, 1], [ y1r,  y0r ], /device
FOR i=0, n2-1 DO BEGIN
  x1b = xbr[i+1]
  x0b = 0.5*( xbr[i] + xbr[i+1] )
  xyouts, x0b, 0.5*(y1r+y0r) - ftsz/2, labs2[i], /device, align=0.5, charsize=2
  plots, x1b*[1, 1], [ y0r, y1r ], /device
ENDFOR

x0b = 0.5*( xbr[0] + xbr[n2] )
xyouts, x0b, plxy - pad/2, 'CLICK MAGNETOGRAM ABOVE TO CYCLE OVERLAY', /device, align=0.5, charsize=2


if (dbl_sw eq 1) then begin

  ;Defining (bottom right) buttons
  labs3 = ['< x10', '<', 'JUMP TO','>',  '> x10']
  n3 = n_elements( labs3 )
  
  xba = 4.0*pad + pls + float(findgen(n3+1))/float(n3)*pls
  
  y0a = pad/4
  y1a = plxy - pad*3/4
  
  ;Plotting revision (bottom left) buttons
  plots, [xba[0], xba[n3]], y0a*[1, 1], /device
  plots, [xba[0], xba[n3]], y1a*[1, 1], /device
  
  plots, xba[0]*[1, 1], [ y1a,  y0a ], /device
  FOR i=0, n3-1 DO BEGIN
    x1b = xba[i+1]
    x0b = 0.5*( xba[i] + xba[i+1] )
    xyouts, x0b, 0.5*(y1a+y0a) - ftsz/2, labs3[i], /device, align=0.5, charsize=2
    plots, x1b*[1, 1], [ y0a, y1a ], /device
  ENDFOR
  
  x0b = 0.5*( xba[0] + xba[n3] )
  xyouts, x0b, plxy - pad/2, 'CLICK MAGNETOGRAM ABOVE TO CYCLE OVERLAY', /device, align=0.5, charsize=2


  ;Defining simultaneous progression
  labs = ['BOTH >', '< BOTH']
  
  n4 = n_elements( labs )
  
  x0d = (1.0 + dbl_sw)*pls+(2.0*(1.0 + dbl_sw) + 0.5)*pad
  x1d = (1.0 + dbl_sw)*pls+(2.0*(1.0 + dbl_sw) - 0.5)*pad+plxx
  ybd = 0.25*pad + float(findgen(n4+1))/float(n4)*(plxy - 0.25*pad)
  
  
  ;Plotting control (right) buttons
  plots, x0d*[1, 1],  [ybd[0], ybd[n4]], /device
  plots, x1d*[1, 1],  [ybd[0], ybd[n4]], /device
  
  plots, [ x1d,  x0d ], ybd[0]*[1, 1], /device
  FOR i=0, n4-1 DO BEGIN
    y1 = ybd[i+1]
    y0 = 0.5*( ybd[i] + ybd[i+1] ) - ftsz/2
    xyouts, 0.5*(x1d+x0d), y0, labs[i], /device, align=0.5, charsize=2
    plots, [ x1d,  x0d ], y1*[1, 1], /device
  ENDFOR



endif


;Identificaiton of buttons
cursor, u, v, 3, /device

;Switch for outside
out_sw = 1

;Cursor inside the control (right) buttons
if ( ( u ge x0c ) and ( u le x1c ) and ( v ge ybc[0] ) and ( v le ybc[nl] ) ) then begin
  ii = where( v LE ybc)
  ii = nl - ii + 1
  out_sw = 0
endif

;Cursor inside the reference (bottom left) buttons
if ( ( u ge xbr[0] ) and ( u le xbr[n2] ) and ( v ge y0r ) and ( v le y1r ) ) then begin
  ii = where( u LE xbr)
  ii = ii + nl
  out_sw = 0
endif

;Cursor inside the left overlay (reference)
if ( ( u ge 2.0*pad ) and ( u le pls+2.0*pad ) and ( v ge pad+plxy ) and ( v le pls+pad+plxy ) ) then begin
  ii = 1 + nl + n2
  out_sw = 0
endif

if (dbl_sw eq 1) then begin

  ;Cursor inside the operations (bottom right) buttons
  if ( ( u ge xba[0] ) and ( u le xba[n3] ) and ( v ge y0a ) and ( v le y1a ) ) then begin
    ii = where( u LE xba)
    ii = ii + nl + n2 + 1
    out_sw = 0
  endif
  
  ;Cursor inside the right overlay (control)
  if ( ( u ge pls+4.0*pad ) and ( u le 2.0*pls+4.0*pad ) and ( v ge pad+plxy ) and ( v le pls+pad+plxy ) ) then begin
    ii = 2 + nl + n2 + n3
    out_sw = 0
  endif
  
  ;Cursor inside the control (right) buttons
  if ( ( u ge x0d ) and ( u le x1d ) and ( v ge ybd[0] ) and ( v le ybd[n4] ) ) then begin
    ii = where( v LE ybd)
    ii = n4 - ii + 3 + nl + n2 + n3
    out_sw = 0
  endif
  

endif

;Cursor inside synchronize button
;if ( ( u ge x0c ) and ( u le x1c ) and ( v ge y0r ) and ( v le y1r ) ) then begin
;  ii = 3 + nl + n2 + n3
;  out_sw = 0
;endif


;Cursor outside buttons area
if  (out_sw eq 1) then begin
  ii = 0
endif

;print, ii

return, ii[0]
END


;-----------------------------------------------------------------------------------------------------------



;EDITS
PRO amj_run_comp, file1=file1, file2 = file2, dbl = dbl, instr

dbl_sw = 0; Switch triggering the double overlay

iskip = 10; number of days to skip

ds_sw = 1; Switch AR overlay

max_sat = 1000;
; constants to change max_sat
incrmnt = 100
max_max_sat = 2000
min_max_sat = 0

if instr eq 3 then begin 

ARs = {ar, mdi_i: 0L, date: '', labl: 0, clr: 0, indxp: '', indxn: '', fluxp: !values.f_nan, fluxn: !values.f_nan, areap:!values.f_nan, arean:!values.f_nan, $
          fcn_ltp: !values.f_nan, fcn_lnp: !values.f_nan, dcenp: !values.f_nan, $ 
          fcn_ltn: !values.f_nan, fcn_lnn: !values.f_nan, dcenn: !values.f_nan, $
          fcenxpp: !values.f_nan, fcenypp: !values.f_nan, dcenpp: !values.f_nan, $
          fcenxpn: !values.f_nan, fcenypn: !values.f_nan, dcenpn: !values.f_nan, $               
          dis: !values.f_nan, tilt: !values.f_nan, lp: !values.f_nan}
endif else begin 

ARs = {ar, mdi_i: 0L, date: '', labl: 0, clr: 0, indxp: '', indxn: '', fluxp: !values.f_nan, fluxn: !values.f_nan, areap:!values.f_nan, arean:!values.f_nan, $
          fcn_ltp: !values.f_nan, fcn_lnp: !values.f_nan, dcenp: !values.f_nan, $ 
          fcn_ltn: !values.f_nan, fcn_lnn: !values.f_nan, dcenn: !values.f_nan, $
          fcenxpp: !values.f_nan, fcenypp: !values.f_nan, dcenpp: !values.f_nan, $
          fcenxpn: !values.f_nan, fcenypn: !values.f_nan, dcenpn: !values.f_nan, $               
          dis: !values.f_nan, tilt: !values.f_nan, lp: !values.f_nan, dm: !values.f_nan, qm: !values.f_nan} 
endelse  

if instr eq 1 then date0 = '1981-06-08';  KPVT 512
if instr eq 2 then date0 = '1992-04-21';  KPVT SPMG
if instr eq 3 then date0 = '1997-01-01';  MDI
if instr eq 4 then date0 = '2010-11-21';  HMI

;First reference day for keeping track time easily
if instr eq 1 then DayOff = julday(1,1,1970); KPVT 512
if instr eq 2 then DayOff = julday(1,1,1970); KPVT SPMG
if instr eq 3 then DayOff = julday(1,1,1993); MDI
if instr eq 4 then DayOff = julday(1,1,2009); HMI

if keyword_set(file1) then begin
  fn1 = file1 
endif 

if (keyword_set(file2) and not keyword_set(file1))then begin
  fn1 = file2 
endif

if (keyword_set(file1) and keyword_set(file2) )then begin
  fn1 = file1
  
  dbl_sw = 1
  
  restore, file2
  PRs2 = PRs
  NRs2 = NRs
  ARs2 = ARs
  mdi_i2 = min(PRs[where(PRs.mdi_i ne 0)].mdi_i)
endif

if ( not keyword_set(file1) and not keyword_set(file2) )then begin  
  fn1 = 'Ar_id_Save0.sav'
endif


restore, fn1
mdi_il = min(PRs[where(PRs.mdi_i ne 0)].mdi_i)

if (dbl_sw eq 1) then begin
  mdi_il = max(mdi_il,mdi_ir)
  mdi_ir = mdi_il
endif

if keyword_set(dbl) then begin

  PRs2 = PRs
  NRs2 = NRs
  ARs2 = ARs
  
  dbl_sw = 1
  mdi_ir = mdi_il

endif

REPEAT BEGIN

  tmp_sw = 0;
  repeat begin
    ;Reading files

    ;Right Magnetogram
    if (dbl_sw eq 1) then begin 
      ;dater = mdi_datestr( string( mdi_ir ), /inverse )
      caldat, mdi_ir + DayOff, Month, Day, Year
      dater = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)
      print, dater
      mgr = amj_file_read( dater, hdrr, instr)
      sr = size(mgr)
      mdi_ir = mdi_ir + 1
      
      if (sr[2] eq 8) then tmp_sw = 1
      
    endif else begin
      tmp_sw = 1
    endelse
    
    ;Left Magnetogram
    ;datel = mdi_datestr( string( mdi_il ), /inverse )
    caldat, mdi_il + DayOff, Month, Day, Year
    datel = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2) 
    mgl = amj_file_read( datel, hdrl, instr)
    sl = size(mgl)
    mdi_il = mdi_il + 1
    
  endrep until ( (tmp_sw eq 1) and (sl[2] eq 8) )
  
  if (dbl_sw eq 1) then mdi_ir = mdi_ir - 1
  mdi_il = mdi_il - 1  
  
  ;Definition of window parameters and window initialization
  pls = 480;
  pad = 40;
  plxy = 100;
  plxx = 200;
  
  sz=size(mgl.img)
  d_zoom = pls/float(sz[1])
  
  window,0,xsize=(1.0 + dbl_sw)*pls+2.0*(1.0 + dbl_sw)*pad+plxx,ysize=pls+2.0*pad+plxy
  
  ;Left Magnetogram
  ;Mixed Overlay
  if (ds_sw eq 1) then begin
    amj_mgplot, mgl.img, mdi_il, instr, PRs = PRs, NRs = NRs, ARs = ARs, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [2.0*pad, pad+plxy, pls+2.0*pad, pls+pad+plxy], title = datel, xrange = [min(mgl.x), max(mgl.x)] , yrange = [min(mgl.y), max(mgl.y)]    
  endif
  ;PNR overlay
  if (ds_sw eq 2) then begin
    amj_mgplot, mgl.img, mdi_il, instr, PRs = PRs, NRs = NRs, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [2.0*pad, pad+plxy, pls+2.0*pad, pls+pad+plxy], title = datel, xrange = [min(mgl.x), max(mgl.x)] , yrange = [min(mgl.y), max(mgl.y)]        
  endif
  ;No overlay
  if (ds_sw eq 3) then begin
    amj_mgplot, mgl.img, mdi_il, instr, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [2.0*pad, pad+plxy, pls+2.0*pad, pls+pad+plxy], title = datel, xrange = [min(mgl.x), max(mgl.x)] , yrange = [min(mgl.y), max(mgl.y)]        
  endif
  
   
  ;Right Magnetogram

  if (dbl_sw eq 1) then begin
  
    ;Mixed Overlay
    if (ds_sw eq 1) then begin
      amj_mgplot, mgr.img, mdi_ir, instr, PRs = PRs2, NRs = NRs2, ARs = ARs2, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [pls+4.0*pad, pad+plxy, 2.0*pls+4.0*pad, pls+pad+plxy], title = dater, xrange = [min(mgr.x), max(mgr.x)] , yrange = [min(mgr.y), max(mgr.y)]    
    endif
    ;PNR overlay
    if (ds_sw eq 2) then begin
      amj_mgplot, mgr.img, mdi_ir, instr, PRs = PRs2, NRs = NRs2, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [pls+4.0*pad, pad+plxy, 2.0*pls+4.0*pad, pls+pad+plxy], title = dater, xrange = [min(mgr.x), max(mgr.x)] , yrange = [min(mgr.y), max(mgr.y)]    
    endif
    ;No overlay
    if (ds_sw eq 3) then begin
      amj_mgplot, mgr.img, mdi_ir, instr, d_xsize = pls, d_ysize = pls, max = max_sat, pos = [pls+4.0*pad, pad+plxy, 2.0*pls+4.0*pad, pls+pad+plxy], title = dater, xrange = [min(mgr.x), max(mgr.x)] , yrange = [min(mgr.y), max(mgr.y)]    
    endif
    
  endif
                       
  stat = 1
  REPEAT BEGIN 
    redraw = 0
    print, 'Please make a choice'
;    tmpcs = 10
;    CASE (tmpcs) of ;button_choice(pls, pad, plxy, plxx) OF
    CASE button_choice(pls, pad, plxy, plxx, dbl_sw) OF
      0: BEGIN
            print, 'Please click inside one of the buttons'
          
         END
      1: BEGIN
            print, 'Quit'
            stat = 0
          
         END
      2: BEGIN
            print, 'Contrast Up'
            max_sat = max_sat - incrmnt
            if (max_sat le min_max_sat) then begin
              max_sat = min_max_sat
              print, 'Maximum contrast reached'
            endif
            print, 'max_sat = ', max_sat
            redraw = 1
  		   END
      3: BEGIN
            print, 'Contrast Down'
            max_sat = max_sat + incrmnt
            if (max_sat ge max_max_sat) then begin
              max_sat = max_max_sat 
             	print, 'Minimum contrast reached'
            endif
            print, 'max_sat = ', max_sat
            redraw = 1
         END
      4: BEGIN
            print, 'Synchronize'
            mdi_il = mdi_ir
            redraw = 1                    
         END     
      5: BEGIN
            print, 'Print Screen'
            tmp_fl = 'prnt_' + strtrim(string(mdi_il,'(I)'),2) + '.png'
            write_png, tmp_fl, TVRD(/true)
         END	  
      6: BEGIN
            print, '<< Left'
            mdi_il = mdi_il - iskip
            redraw = 1
         END
      7: BEGIN
            print, '< Left'
            mdi_il = mdi_il - 1
            redraw = 1
         END
  	  8: BEGIN
            print, 'Jump to Left'
      			date = ''
      			read, 'Enter date (e.g. 2001-08-31): ', date
      			mdi_tmp = long( mdi_datestr( date ) )
      			date_tmp = mdi_datestr( string( mdi_tmp ), /inverse ) 
      			mg_tmp = amj_file_read( date_tmp, hdrr, instr )
      			
      			tmpsize = size(mg_tmp)
      			if (tmpsize[2] ne 8) then begin
      				print, 'No available MDI magnetogram.'
      			endif else begin
      				mdi_il = mdi_tmp
      				redraw = 1
      			endelse  			
  		   END
      9: BEGIN
            print, '> Left'
            mdi_il = mdi_il + 1
            redraw = 1
         END
     10: BEGIN
            print, '>> Left'
            mdi_il = mdi_il + iskip
            redraw = 1          
         END
     11: BEGIN
            print, 'Left Overlay'
            ds_sw = ds_sw + 1
            if (ds_sw gt 3) then ds_sw = 1
            redraw = 1                    
         END
     12: BEGIN
            print, '<< Right'
            mdi_ir = mdi_ir - iskip
            redraw = 1
         END
     13: BEGIN
            print, '< Right'
            mdi_ir = mdi_ir - 1
            redraw = 1
         END
     14: BEGIN
            print, 'Jump to Right'
            date = ''
            read, 'Enter date (e.g. 2001-08-31): ', date
            mdi_tmp = long( mdi_datestr( date ) )
            date_tmp = mdi_datestr( string( mdi_tmp ), /inverse ) 
            mg_tmp = amj_file_read( date_tmp, hdrr, instr )
            
            tmpsize = size(mg_tmp)
            if (tmpsize[2] ne 8) then begin
              print, 'No available MDI magnetogram.'
            endif else begin
              mdi_ir = mdi_tmp
              redraw = 1
            endelse       
         END
     15: BEGIN
            print, '> Right'
            mdi_ir = mdi_ir + 1
            redraw = 1
         END
     16: BEGIN
            print, '>> Right'
            mdi_ir = mdi_ir + iskip
            redraw = 1          
         END
     17: BEGIN
            print, 'Right Overlay'
            ds_sw = ds_sw + 1
            if (ds_sw gt 3) then ds_sw = 1
            redraw = 1                    
         END
     18: BEGIN
            print, '< Both'
            mdi_ir = mdi_ir - 1
            mdi_il = mdi_il - 1
            redraw = 1
         END
     19: BEGIN
            print, 'Both >'
            mdi_ir = mdi_ir + 1
            mdi_il = mdi_il + 1
            redraw = 1
         END
    ENDCASE
	
  ENDREP UNTIL( redraw or ( stat eq 0 ) )
ENDREP UNTIL ( stat EQ 0 )

return
END


;----------------------------------------------------------------------------------------------------------

 
