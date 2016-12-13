pro mksrc3,file,l,b,gsed_file=gsed_file,map_ext=map_ext,sed_ext=sed_ext


sed=mrdfits(gsed_file,sed_ext)

spec_gsed=sed.flux		
ene_gsed=sed.energy

lene_gsed=alog10(ene_gsed/1d12)

lene=([findgen(300)/100.]-1.)
ene=10.d0^lene                                   ;TeV

spec=interpol(spec_gsed,lene_gsed,lene)/(1.6)/ene^2.     ;  ph/cm sec TeV

lemin=findgen(30)*0.1-1. 
lemax=lemin+0.1
flux=fltarr(n_elements(lemin))

for i=0,n_elements(lemin)-1 do begin
       
      w=where(lene ge lemin(i) and lene lt lemax(i))
      flux(i)=meag(spec(w))*(10.^lemax(i)-10.^lemin(i))

endfor
		
		
		
;;; >>>>>> Scrivi su file (if:diffuse/else:poit-like) <<<<<<

if keyword_set(map_ext) then begin
	im = mrdfits(gsed_file, map_ext, im_header)
	
	map_size=size(im)
	sky=fltarr(map_size[1],map_size[2],n_elements(flux))
 
        for i=0,n_elements(flux)-1 do begin
		sky[*,*,i]=im*flux[i] / total(im)
	endfor

	MKHDR,he1,sky,/image

	com=''
	
	fxaddpar,he1,'CRVAL1',fxpar(im_header,'CRVAL1',comment=com),com 
	fxaddpar,he1,'CRPIX1',fxpar(im_header,'CRPIX1',comment=com),com 
	fxaddpar,he1,'CDELT1',fxpar(im_header,'CDELT1',comment=com),com 
	fxaddpar,he1,'CTYPE1',fxpar(im_header,'CTYPE1'),' Axis type for dim 1'
	fxaddpar,he1,'CUNIT1','deg     ',''
	fxaddpar,he1,'LTM1_1',fxpar(im_header,'LTM1_1',comment=com),com

	fxaddpar,he1,'CRVAL2',fxpar(im_header,'CRVAL2',comment=com),com 
	fxaddpar,he1,'CRPIX2',fxpar(im_header,'CRPIX2',comment=com),com 
	fxaddpar,he1,'CDELT2',fxpar(im_header,'CDELT2',comment=com),com
	fxaddpar,he1,'CTYPE2',fxpar(im_header,'CTYPE2'),' Axis type for dim 2'
	fxaddpar,he1,'CUNIT2','deg     ',''
	fxaddpar,he1,'LTM2_2',fxpar(im_header,'LTM2_2',comment=com),com

	fxaddpar,he1,'EQUINOX',2000.0,' Equinox of Ref. Coord.'
	fxaddpar,he1,'RADECSYS','FK5',' WCS for this file'

	fxaddpar,he1,'CRVAL3',1,''
	fxaddpar,he1,'CRPIX3',0,''
	fxaddpar,he1,'CDELT3',1,''

	fxaddpar,he1,'PIXCENT','T','' 
	fxaddpar,he1,'BUNIT','INTENSITY',''
	fxaddpar,he1,'BSCALE',fxpar(im_header,'BSCALE',comment=com),com 
	fxaddpar,he1,'BZERO', fxpar(im_header,'BZERO',comment=com),com
	fxaddpar,he1,'PRIMTYPE','RECT_MAP',''
	fxaddpar,he1,'EXTNAME','FLUX'


	writefits,file,sky,he1,/append

endif else begin
	sky=fltarr(n_elements(l),n_elements(b),n_elements(flux))
	sky(*)=flux

	if n_elements(l) eq 1 then begin
	bin=1.
	endif else bin=l(1)-l(0)

	MKHDR,he1,sky,/image

	fxaddpar,he1,'CRVAL1',l(0),''
	fxaddpar,he1,'CRPIX1',1,''
	fxaddpar,he1,'CDELT1',-bin,''
	fxaddpar,he1,'CTYPE1','GLON-CAR',''
	;fxaddpar,he1,'CTYPE1','GLON',''
	fxaddpar,he1,'CUNIT1','deg     ',''

	fxaddpar,he1,'CRVAL2',b(0),''
	fxaddpar,he1,'CRPIX2',1,''
	fxaddpar,he1,'CDELT2',bin,''
	fxaddpar,he1,'CTYPE2','GLAT-CAR',''
	;fxaddpar,he1,'CTYPE2','GLAT',''
	fxaddpar,he1,'CUNIT2','deg     ',''

	fxaddpar,he1,'CRVAL3',1,''
	fxaddpar,he1,'CRPIX3',0,''
	fxaddpar,he1,'CDELT3',1,''

	giu:fxaddpar,he1,'PIXCENT','T',''
	fxaddpar,he1,'BUNIT','INTENSITY',''
	fxaddpar,he1,'BSCALE',1.,''
	fxaddpar,he1,'BZERO',0.,''
	fxaddpar,he1,'PRIMTYPE','RECT_MAP',''
	fxaddpar,he1,'EXTNAME','FLUX'

	writefits,file,sky,he1,/append

      endelse

;;;;  Energy Bounds

FXBHMAKE,head,n_elements(lemin),'EBOUNDS'

fxbaddcol,1,head,float(1),'L_EMIN'
fxbaddcol,2,head,float(1),'L_EMAX'
fxbcreate,unit,file,head

FOR I=0,n_elements(lemin)-1. do begin
  
  fxbwrite,unit,lemin(i),1,I+1
  fxbwrite,unit,lemax(i),2,I+1
  
 ENDFOR

fxbfinish,unit

print,'file: "',file,'" >>> Done!'

end

