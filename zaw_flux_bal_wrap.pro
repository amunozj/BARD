PRO zaw_flux_bal_wrap, instr

SET_PLOT, 'Z'
thisDevice = !D.Name
dir = './Saves/*.sav'
spawn, 'ls ' + dir, result

for i = 0, n_elements(result) - 1 do begin 
    print, result[i]
    zaw_save_flux_bal, instr, result[i], O_fldr = './Balanced/'
endfor
    
    
return
END