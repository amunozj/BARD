PRO amj_pnr_pre_frag_wrap

SET_PLOT, 'Z'
thisDevice = !D.Name

dir = './Pre_dt_brk/tmp/*.sav'
spawn, 'ls ' + dir, result

for i = 0, n_elements(result)-1 do begin  

    print, result[i]
    amj_pnr_pre_frag, result[i], O_fldr='./Pre_dt_brk/tmp/'    

endfor
    
    
return
END


;----------------------------------------------------------------------------------------------------------
