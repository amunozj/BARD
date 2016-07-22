;-------------------------------------------------------------------------------------------------------------------------------
; NAME:       zaw_save_flux_bal
;
; PURPOSE:	Inputs a save file and operates on each day, balancing flux with zaw_ar_flux_bal
;
; CALLING SEQUENCE: zaw_ar_flx_balance, instr, input_file
;
; INPUTS:		instr 				Instrument of save file for default settings
;				input_file			Location of file to be processed
;
; OPTIONAL INPUTS: 0_fldr			Location to output the save file. Must be in format './directory/'
;
;-------------------------------------------------------------------------------------------------------------------------------
PRO zaw_save_flux_bal, instr, input_file, O_fldr = O_fldr

;Does not plot onscreen
SET_PLOT, 'Z'
thisDevice = !D.Name


;Initializing run and variables
;-------------------------------------------------------------------------------------

;Region detection constants

;MDI
;seg_const={ker_th:500.0, ker_th2:275.0, ar_th:325.0, ar_th2:50.0, eros_size:9.0, dila_size:18.0, npssu: 4, dis_lim:3.0, ovr_lim: 0.3, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}

;KPVT
seg_const={ker_th:400.0, ker_th2:200.0, ar_th:150.0, ar_th2:50.0, eros_size:10.0, dila_size:20.0, npssu:2, dis_lim:2.0, ovr_lim: 0.2, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}


;Instrument Specifications
if instr eq 1 then DayOff = julday(1,1,1970);	KPVT 512
if instr eq 2 then DayOff = julday(1,1,1990);	KPVT SPMG
if instr eq 3 then DayOff = julday(1,1,1993);	MDI
if instr eq 4 then DayOff = julday(1,1,2009);	HMI

;Active region detection constants
ar_cnst={dis_lim1:5.0, dis_lim2:4.0, exp_f: 1.0,  exp_d: 4.0, exp_s: 1.0, mlth: 40.0, mxB: 180.0, MxFlxim:3.0, Imb_tol: 0.10, Imb_it: 10, lim_lon: -90.0, k_sig:15.0, npr: 5, nmgnt: 5 , vld_thr: 0.69, valid_range:[-20000.,20000.]}

;KPVT 512
if instr eq 1 then begin
	seg_const={ker_th:400.0, ker_th2:200.0, ar_th:150.0, ar_th2:50.0, eros_size:10.0, dila_size:20.0, npssu:2, dis_lim:2.0, ovr_lim: 0.2, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}
endif

;KPVT SMPG
if instr eq 2 then begin
	seg_const={ker_th:400.0, ker_th2:200.0, ar_th:150.0, ar_th2:50.0, eros_size:10.0, dila_size:20.0, npssu:2, dis_lim:2.0, ovr_lim: 0.2, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}
endif

;MDI
if instr eq 3 then begin
	seg_const={ker_th:500.0, ker_th2:275.0, ar_th:325.0, ar_th2:50.0, eros_size:9.0, dila_size:18.0, npssu: 3, dis_lim:3.0, ovr_lim: 0.3, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}
endif

;HMI
if instr eq 4 then begin
	seg_const={ker_th:200.0, ker_th2:200.0, ar_th:150.0, ar_th2:50.0, eros_size:10.0, dila_size:30.0, npssu:1, dis_lim:2.0, ovr_lim: 0.0, qr_th:60.0, mxB: 0.0, k_sig:15.0, ar_grow_sig:20., valid_range:[-20000.,20000.], deg_lim:70.0}
endif

restore, input_file
output_file = input_file

mdi_ir  = min([min(PRs[where(PRs.mdi_i ne 0)].mdi_i),min(NRs[where(NRs.mdi_i ne 0)].mdi_i)])
mdi_if  = max([max(PRs[where(PRs.mdi_i ne 0)].mdi_i),max(NRs[where(NRs.mdi_i ne 0)].mdi_i)])

caldat, mdi_ir+DayOff, Month, Day, Year
start_date = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)

caldat, mdi_if+DayOff, Month, Day, Year
end_date = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)

while (mdi_ir le mdi_if) do begin

;Reading files
repeat begin
    caldat, mdi_ir+DayOff, Month, Day, Year
	dater = strtrim(string(Year),2)+'-'+strtrim(string(Month,format='(I02)'),2)+'-'+strtrim(string(Day,format='(I02)'),2)
	print, 'Balancing ' + dater
	mgr = amj_file_read( dater, hdrr, instr )
	          
	sr = size(mgr)

    mdi_ir = mdi_ir + 1  
endrep until ( (sr[2] eq 8) or (mdi_ir gt mdi_if) )
mdi_ir = mdi_ir - 1 		;to reset day back to normal

ar_in = where((ARs.mdi_i eq mdi_ir), n_ars)
if (n_ars ne 0) then begin
	amj_coord, mgr.img, hdrr, CRD, instr, seg_const=seg_const 
	zaw_ar_flux_bal, CRD, mdi_i, ARs, ar_in, n_ars, dater, ar_cnst = ar_cnst, seg_const = seg_const
endif
mdi_ir = mdi_ir + 1

endwhile 

if keyword_set(O_fldr) then begin
      output_fn = O_fldr + 'Ar_' + start_date + '_to_' + end_date + '_flux_bal.sav'
endif else begin
	output_fn = 'Ar_' + start_date + '_to_' + end_date + '_flux_bal.sav'
endelse

print, "Saving " + output_fn
SAVE, ARs, PRs, NRs, PRsFr, NRsFr, mdi_il, mdi_ir, lbl, prepnr_sw, mdi_ir_vis, buff_sw, instr, FILENAME = output_fn
print, "Saved"
END