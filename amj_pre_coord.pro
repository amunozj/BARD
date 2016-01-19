PRO amj_pre_coord, instr, start_date=start_date, end_date=end_date

SET_PLOT, 'Z'
thisDevice = !D.Name

; Setting file structure-----------------------------------------------------------

dir0= '/disk/data/munoz/BARD/precoord/'

; KPVT 512
if instr eq 1 then begin
  dir_str = 'KPVT'
  seg_const={ker_th:400.0, ker_th2:200.0, ar_th:150.0, ar_th2:50.0, eros_size:10.0, dila_size:20.0, npssu:2, dis_lim:2.0, ovr_lim: 0.2, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}
endif

;KPVT SMPG
if instr eq 2 then begin
  dir_str = 'SMPG'
  seg_const={ker_th:400.0, ker_th2:200.0, ar_th:150.0, ar_th2:50.0, eros_size:10.0, dila_size:20.0, npssu:2, dis_lim:2.0, ovr_lim: 0.2, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}
endif 
  
;MDI
if instr eq 3 then begin
  dir_str = 'MDI'
  seg_const={ker_th:500.0, ker_th2:275.0, ar_th:325.0, ar_th2:50.0, eros_size:9.0, dila_size:18.0, npssu: 3, dis_lim:3.0, ovr_lim: 0.3, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}
endif


;HMI
if instr eq 4 then begin
  dir_str = 'HMI'
  seg_const={ker_th:200.0, ker_th2:200.0, ar_th:150.0, ar_th2:50.0, eros_size:10.0, dila_size:30.0, npssu:1, dis_lim:2.0, ovr_lim: 0.0, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}
endif

dir = dir0 + dir_str + '/'


;First reference day for keeping track time easily
if instr eq 1 then DayOff = julday(1,1,1970); KPVT 512
if instr eq 2 then DayOff = julday(1,1,1970); KPVT SPMG
if instr eq 3 then DayOff = julday(1,1,1993); MDI
if instr eq 4 then DayOff = julday(1,1,2009); HMI

mdi_ir  = julday(strmid(start_date,5,2),strmid(start_date,8,2),strmid(start_date,0,4))-DayOff
mdi_if  = julday(strmid(end_date,5,2),strmid(end_date,8,2),strmid(end_date,0,4))-DayOff
 
while (mdi_ir le mdi_if) do begin

  repeat begin
    ;Reading files
    caldat, mdi_ir+DayOff, Month, Day, Year
    dater = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)
    mgr = amj_file_read( dater, hdrr, instr )    
            
    sr = size(mgr)

    mdi_ir = mdi_ir + 1
    print, mdi_ir
    
  endrep until ( (sr[2] eq 8) or (mdi_ir gt mdi_if) )  
  mdi_ir = mdi_ir - 1

  amj_coord, mgr.img, hdrr, CRD, instr, seg_const=seg_const

  SAVE, CRD, FILENAME = dir + 'Pre_Coord_' +  dir_str + '_' + strtrim(string(Year),2) + strtrim(string(Month,format='(I02)'),2) + strtrim(string(Day,format='(I02)'),2) + '.sav'
  
  mdi_ir = mdi_ir + 1    
endwhile

return
END
