PRO amj_pre_coord_read, instr, mdi_ir, CRD

print, 'Reading pre-calculated coordinate information'

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
if instr eq 2 then DayOff = julday(1,1,1990); KPVT SPMG
if instr eq 3 then DayOff = julday(1,1,1993); MDI
if instr eq 4 then DayOff = julday(1,1,2009); HMI

;Reading file
caldat, mdi_ir+DayOff, Month, Day, Year

restore, dir + 'Pre_Coord_' +  dir_str + '_' + strtrim(string(Year),2) + strtrim(string(Month,format='(I02)'),2) + strtrim(string(Day,format='(I02)'),2) + '.sav'
  
return
END
