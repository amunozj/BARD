pro mreadfits_tilecomp, tfiles, index, data, $
   llxp, llyp, nxp, nyp, _extra=_extra, $
   status=status, old_paradigm=old_paradigm, $
   use_index=use_index, fovpix=fovpix, outdir=outdir, parent_out=parent_out, $
   funcomp=funcomp, silent=silent, chrontree=chrontree, $ 
   comp_delete=comp_delete, uncomp_delete=uncomp_delete, $
   fnames_uncomp=fnames_uncomp, only_uncompress=only_uncompress, $
   noshell=noshell, debug=debug, time_tag=time_tag, $
   use_shared_lib=use_shared_lib,shared_lib_path=shared_lib_path
;
;   Purpose: read tile compressed fits files (used by jsoc/sdo for example)
;
;   Input Parameters:
;      tfiles - file(s) to read
;      index - if 'index' input is structure
;      xllp, yllp, nxp, nyp - optional subfield, in pixels
;                             lower left xllp & yllp may be vectors 1:1 files
;
;   Output Parameters:
;     
;      data - 2D/3D image array 
;
;   Keyword Paramters:
;      use_index - if 'index' is input & structure and /USE_INDEX, then...
;      comp_delete - if set, remove input/compresssed versions
;      uncomp_delete if set, remove output/Uncompressed (only if DATA present)
;      silent - if set, inhibit some status/diagnostics
;      fovpix - optional subfield in pixels (keyword synonym for postional) [llx,lly,nx,ny]
;      parent_out - path for uncompressed files (default=get_temp_dir() )
;      outdir - synonym for PARENT_OUT
;      _extra - keyword inherit -> mreadfits_shm or mreadfits
;     fnames_uncomp (output) - full paths to uncompressed files
;     only_uncompress - no data read, just uncompress and return FNAMES_UNCOMP
;     chrontree - if set, ucompressed files -> "standard" chronologically
;                 organized tree like: <parent>/yyyy/mm/dd/Hhhhh/<files>
;     noshell (switch) - if set, attempt to avoid shell w/imcopy (full image only)
;     use_shared_lib (switch) - use call_external interface to cfitsio shared library (instead of
;                               spawning imcopy), note this is OS dependent and possibly IDL version
;                               dependent, and has the potential to crash IDL if the shared library
;                               is not properly compiled
;     shared_lib_path - optional keyword specifying path of shared library (with trailing slash) for
;                       fitsio.so, default is the same path as for imcopy
;
;   History:
;      4-mar-2010 - S.L.Freeland 
;      8-mar-2010 - S.L.Freeland - post imcopy -> distrib - proto working
;      8-apr-2010 - S.L.Freeland - subfield support via xllp,yllp,nxp,nyp (in pixels)[
;     22-jul-2010 - S.L.Freeland - avoid an issue 
;     10-aug-2010 - S.L.Freeland - change default PARENT_OUT to get_temp_dir()
;     10-sep-2010 - S.L.Freeland - assure that uncompressed output filenams are uniq
;              			   add OUTDIR synonym for PARENT_OUT
;     18-oct-2010 - S.L.Freeland - made /SILENT silenter
;     24-jan-2011 - S.L.Freeland - add /NOSHELL keyword & tweak (full image only)
;     31-mar-2011 - S.L.Freeland - allow vectorized FOV - removed
;                   /NOSHELL restriction for FOV (but caveat emptor)
;     17-may-2011 - S.L.Freelnd - fix typo-bug with NY
;                   write_access.pro -> file_test(/write) intrinsic
;     25-aug-2011 - avoid collision between /USE_INDEX and new tile paradigm
;     18-Apr-2012 - Zarro (ADNET), pass _extra to spawn
;     27-sep-2012 - M.L.DeRosa - added shared library option for uncompression
;      2-oct-2012 - S.L.Freeland - respect /silent flag a bit more
;     23-oct-2012 - S.L.Freeland - merge {freeland/zarro/derosa} changes
;                  use ssw_bin_path.pro for imcopy & fitsio.so OS/ARCH lookup
;      5-apr-2013 - S.L.Freeland - "missing" tag protection (IRIS)
;     18-jul-2013 - S.L.Freeland - call mreadfits.pro for mixed nx/ny (caveat emptor)
;
;   Restrictions:
;      In theory, this will evolve to a "gen" routine - but for now,
;      a little SDO AIA/HMI centric.
;      as of "now", only a subset of subfield fits keywords reflect extraction (todo in wrapper?)
;
;   Method:
;     OS/ARCH dependent cfitsio library (imcopy for today)
;
;   Warning:
;     Potential disk killer - pay attention to PARENT_OUT (aka OUTDIR) and /UNCOMP_DELETE options
;     /NOSHELL will usually speed this up -BUT- - not all os/arch/shells
;        tested - try it, if it works, good, if not, leave it off...
;
;-
; if n_elements(index) ne n_elements(tfiles) then $
old_paradigm=keyword_set(old_paradigm)
old_paradigm=old_paradigm or (keyword_set(use_index) and n_elements(index) eq n_elements(tfiles))
noshell=keyword_set(noshell) ; and n_params() le 3  ; 31-mar - removed FOV


loud=1-keyword_set(silent)
if not old_paradigm then $
   mreadfits_header,tfiles,index,only_tags=only_tags,extension=1

use_shared_lib=keyword_set(use_shared_lib) or data_chk(shared_lib_path,/string)

if use_shared_lib then begin 
   defpath=ssw_bin_path('fitsio.so', found=found,/path_only,/ontology)
   case 1 of 
      keyword_set(shared_lib_path): so_path=shared_lib_path ; user supplied
      else: so_path=defpath ; default path
   endcase
   found=file_exist(concat_dir(so_path,'fitsio.so')) ; verify available for OS/ARC
   if not found then begin 
      box_message,'fitsio/shared object request but not available for this OS/ARCH - using imcopy
      use_shared_lib=0 ; override
   endif
endif

imcopy=ssw_bin_path('imcopy', /ontology, found=imcfound)

ndata=n_params() gt 2 and ~keyword_set(only_uncompress)

case 1 of 
   n_params() lt 2: begin
      box_message,'Need at least TFILES & INDEX input'
      return
   endcase
   else: 
endcase
nfiles=n_elements(tfiles)

if not file_exist(tfiles(0)) then begin 
   box_message,'File not found, returning...'
   return
endif 

htest=headfits(tfiles(0),ext=1,errmsg=errmsg) 
if errmsg(0) ne '' then begin 
   box_message,['Error reading file',errmsg]
   return
endif

allhead=strarr(nfiles * n_elements(htest))  ; header buffer 

shtest=fitshead2struct(htest)

if ndata and n_params() lt 5 then begin ; need DATA output array
   nx=gt_tagval(shtest,/znaxis1,missing=-1)
   ny=gt_tagval(shtest,/znaxis2,missing=-1)
   if nx eq -1 or ny eq -1 then begin
      box_message,'Need data but cant derive NX/NY from extended header...'
      return
   endif
   data=make_array(nx,ny,nfiles,/uint,/nozero)
endif

c2u=str_subset(shtest,'zbitpix,znaxis,znaxis1,znaxis2,blank,bzero,bscale')
tc2u=tag_names(c2u)
ztags=where(strmid(tc2u,0,1) eq 'Z',zcnt)
if zcnt gt 0 then begin
   zt=tc2u(ztags)
   nozt=strmid(zt,1,8)
   for i=0,zcnt-1 do c2u=rep_tag_name(c2u,zt(i),nozt(i))
end

c2u=replicate(c2u,nfiles)
index=join_struct(c2u,index)

; check for subfield request... (31-mar-2011 - recast for 3D)
if n_params() ge 6 then begin 
   if n_elements(llxp) eq nfiles then $
     xll=llxp else xll=replicate(llxp,nfiles) ; vector ok for XLL/YLL
   if n_elements(llyp) eq nfiles then $
     yll=llyp else yll=replicate(llyp,nfiles)
   nx=replicate(nxp(0),nfiles)  ; scalar only for NX/NY
   if n_elements(nyp) eq 0 then ny=nx else ny=replicate(nyp(0),nfiles)

   sfov='['+ $
            strtrim(xll,2) + ':' + strtrim(xll+nx-1,2) + ',' + $
            strtrim(yll,2) + ':' + strtrim(yll+ny-1,2) +       $
         ']'
endif
fov=''

;  branch if you're daring enough to try using the shared library...
if keyword_set(use_shared_lib) then begin

  ;  define data array
  case gt_tagval(shtest,/zbitpix,missing=-1) of 
          8:   IDL_type = 1          ; Byte
         16:   IDL_type = 2          ; Integer*2
         32:   IDL_type = 3          ; Integer*4
         64:   IDL_type = 14         ; Integer*8
        -32:   IDL_type = 4          ; Real*4
        -64:   IDL_type = 5          ; Real*8
       else:   begin
               message,/CON, 'ERROR - Illegal value of BITPIX (= ' +  $
               strtrim(bitpix,2) + ') in FITS header'
               return
               end
  endcase     
  data=0  ;  saves memory in some cases when data has been defined above
  data=make_array(nx(0),ny(0),nfiles,/nozero,type=IDL_type)

  ;  read in data using shared library
  tfilesq=tfiles+(['','[1]'])(n_elements(headfits(tfiles(0),/exten)) gt n_elements(headfits(tfiles(0))))
  if n_elements(sfov) eq nfiles then tfilesq+=sfov
  for i=0,nfiles-1 do begin
    if loud then box_message,'reading '+tfiles(i)
    data(*,*,i)=fitsio_read_image(tfilesq(i),so_path=so_path)
  endfor
    
endif else begin  ;  ... or use the tried-and-true imcopy utility

  tempdir=get_temp_dir()
  pdir=tempdir
  if keyword_set(outdir) then parent_out=outdir(0) ; synonyms
  if data_chk(parent_out,/string) then begin 
     if file_test(parent_out,/write) then pdir=parent_out else $
        box_message,'No write access to PARENT_OUT='+parent_out
  endif

  if loud then box_message,'Uncompressing to> ' + pdir
  if get_logenv('check_time') ne '' then stop, 'time_tag'
  fnames=ssw_jsoc_index2filenames(index,parent_out=pdir, $
     chrontree=chrontree,/mkdir, time_tag=time_tag)
  nnames=n_elements(fnames)
  unss=uniq(fnames)
  ; block added 10-aug-2010 to assure Uniq uncompressed file names
  if n_elements(unss) lt nnames then begin 
     xexts=ssw_strsplit(tfiles,'.',/tail,head=xfnames)
     xfnames=strmids(xfnames,str_lastpos(xfnames,'/')+1,1000)
     ucss=uniq(xfnames,sort(xfnames))
     if n_elements(ucss) eq n_elements(xfnames) then fnames=concat_dir(pdir,xfnames+'_uncomp.'+xexts) else begin 
       box_message,'still not uniq... xfnames'
       suffix='_'+string(sindgen(nnames),format='(i3.3)')
       fnames=concat_dir(pdir,xfnames+suffix+'.'+xexts)
     endelse
  endif

  debug=keyword_set(debug)
  fnames_uncomp=fnames

  ; ~transparent handle of primary vs extension image location
  tfilesq=tfiles+(['','[1]'])(n_elements(headfits(tfiles(0),/exten)) gt n_elements(headfits(tfiles(0))))
  if n_elements(sfov) eq nfiles then begin
     imccmds=imcopy + " '" + tfilesq+sfov + "' " + fnames ; subfield requires quoting ; TODO WinXX?
  endif else imccmds=imcopy + " '" + tfilesq +"' '" + fnames + "'"
  if noshell then imccmds=str_replace(imccmds,"'","") ; remove single quotes around filenames

  if debug then stop,'imccmds'

  for i=0, nfiles-1 do begin
     if loud then box_message,tfiles(i) + ' -> ' + fnames(i)
     ssw_file_delete,fnames(i)
     ; spawn,str2arr(imccmds(i),' '),/noshell 
     if noshell then spawn,str2arr(imccmds(i),' '),/noshell,_extra=_extra else $ ; attempt; no guarantees...
        spawn,imccmds(i),_extra=_extra ; for now at least, must use shell...
     if get_logenv('noshell') ne '' then stop,'noshell check'
     if keyword_set(old_paradigm) then begin  
        mreadfits_shm,fnames(i),index(i),idata,/noscale
        if get_logenv('check') then stop,'fnames(i)'
        temp=index(i)
        mreadfits_fixup,temp,idata
        if required_tags(index,'cdelt1,cdelt2') then begin 
           temp.cdelt1=index(i).cdelt1
           temp.cdelt2=index(i).cdelt2
        endif
        index(i)=temp
        ssw_file_delete,fnames(i)
        mwritefits,index(i),idata,outfile=fnames(i)
     endif
  endfor

  uncomp_delete=keyword_set(uncomp_delete) or required_tags(_extra,'uncomp_delete') ; 
  if ndata then begin 
     znaxis1=gt_tagval(index,/znaxis1,missing=1)
     znaxis2=gt_tagval(index,/znaxis2,missing=1)
     if n_elements(all_vals(znaxis1)) gt 1 or n_elements(all_vals(znaxis2)) gt 1 then begin 
        if loud then box_message,'mixed nx/ny - embedding in largest nx/ny; caveat emptor'
        mreadfits,fnames,index,data,_extra=_extra,/insert
     endif else mreadfits_shm,fnames,index,data,_extra=_extra
  endif else if loud then box_message,'no data read'

  if  uncomp_delete then  begin 
     if loud then box_message,'Removing uncompressed versions on request'
     ssw_file_delete, fnames ; user requests delete
  endif

endelse

return
end



