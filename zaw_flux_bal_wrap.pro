;-------------------------------------------------------------------------------------------------------------------------------
; NAME:      	zaw_flux_bal_wrap
;
; PURPOSE:      Finds .sav files and applies the zaw_save_flux_bal routine
;
; CALLING SEQUENCE: zaw_flux_bal_wrap, instr
;
;  
; INPUTS:			instr						Instrument of save file for default settings
;
; OPTIONAL INPUTS:	dir_in					Location of the directory of where the saves will be located. Must be in format './directory/'
;					dir_out					Location to output the save files. Must be in format './directory/'
;
;-------------------------------------------------------------------------------------------------------------------------------
PRO zaw_flux_bal_wrap, instr, dir_in = dir_in, dir_out = dir_out

SET_PLOT, 'Z'
thisDevice = !D.Name

if keyword_set(dir_in) then begin
	dir = dir_in + '*.sav'
	spawn, 'ls ' + dir, result
endif else begin
	spawn, 'ls *.sav', result
endelse

for i = 0, n_elements(result) - 1 do begin 
    print, result[i]
    if keyword_set(dir_out) then begin
    	zaw_save_flux_bal, instr, result[i], O_fldr = dir_out
    endif else begin
    	zaw_save_flux_bal, instr, result[i]
    endelse
endfor
      
return
END