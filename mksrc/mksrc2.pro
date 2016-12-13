pro mksrc2,file,l,b,emin,emax,z_src=z_src,NEbins=NEbins,spec_type=spec_type,k=k,index1=index1,index2=index2, $
		E0=E0,Ecut=Ecut,Emean=Emean,Sigma=Sigma,spec_file=spec_file,DiffMap=DiffMap,EXTEN_NO=EXTEN_NO

; The script creates MapCube of the point-like/diffuse source
;
; file					: output MapCube
; (l,b)					: galactic coordinates in deg
; (emin,emax)				: minimum/maximum limits of the energy range (TeV)
; z_src         [default=0.00]          : redshift of the source 
; NEbins 	[default=30]		: number of energy bins 
; spec_type 	[default='PowerLaw']	: spectral type ('ConstantValue','PowerLaw','PowerLaw2','ExpCutoff', $
;					   		 'BrokenPowerLaw','LogParabola','Gaussian','FileFunction') 
; k 		[default=1.0e-11]	: flux @ E0 TeV in [ph cm-2 s-1 TeV-1] for all spec_type except $
;								'PowerLaw2' and 'Gaussian' in [ph cm-2 s-1]
; index1 	[default=-2.0]		: spectral index
; index2	[default=-3.0] 		: second spectral index/curvature for spec_type='BrokenPowerLaw'/'LogParabola' 
; E0		[default=1.0 TeV]	: scale energy or break energy for spec_type='BrokenPowerLaw' 
; Ecut		[default=10.0 TeV]	: cut-off energy for spec_type='ExpCutoff'
; Emean		[default=5.0 TeV]	: Mean energy for spec_type='Gaussian'
; Sigma		[default=1.0 TeV]	: Sigma for spec_type='Gaussian'
; spec_file				: file for spec_type='FileFunction', an ASCII file with 2 columns of energy [TeV] and $
;										differential flux [ph cm-2 s-1 TeV-1]
; DiffMap				: diffuse-map file of the extended source
; EXTEN_NO				: diffuse-map extension (number)

if ((n_elements(file) ne 0) && (n_elements(l) ne 0) && (n_elements(b) ne 0) && $
		(n_elements(emin) ne 0) && (n_elements(emax) ne 0)) then begin
	
	if not keyword_set(NEbins) then NEbins=30
        if not keyword_set(z_src)  then z_src=0.0
	if not keyword_set(spec_type) then spec_type='PowerLaw'
	if not keyword_set(index1) then index1=-2.0
	if not keyword_set(index2) then index2=-3.0
	if not keyword_set(k) then k=1.0e-11
	if not keyword_set(E0) then E0=1.0
	if not keyword_set(Ecut) then Ecut=10.0
	if not keyword_set(Emean) then Emean=5.0
	if not keyword_set(Sigma) then Sigma=1.0
	print, file,':',l,b,emin,emax,z_src,NEbins,' ',spec_type
	print,' '

	n=10 ; number of flux point in each energy bin
	log_lemin = alog10(emin)
	log_lemax = alog10(emax)
	e_incr = (log_lemax-log_lemin)/(NEbins*n)
	lene = [log_lemin:log_lemax:e_incr]
	ene = 10.d0^lene
	;;; List of spectral types	
	case 1 of
		spec_type eq 'ConstantValue'	: spec=k*ene/ene
		spec_type eq 'PowerLaw'		: spec=k*(ene/E0)^index1
		spec_type eq 'PowerLaw2'	: spec=k*(index1+1)*ene^index1/(lemax^(index1+1)-lemin^(index1+1))
		spec_type eq 'ExpCutoff'	: spec=k*(ene/E0)^index1*exp(-ene/Ecut)
		spec_type eq 'BrokenPowerLaw'	: begin
							spec=[]
							if E0 ge lemin then spec=[spec,k*(ene[where(ene le E0)]/E0)^index1]
							if E0 le lemax then spec=[spec,k*(ene[where(ene gt E0)]/E0)^index2]
       						  end
		spec_type eq 'LogParabola'	: spec=k*(ene/E0)^(index1+index2*alog(ene/E0))
		spec_type eq 'Gaussian'		: spec=k/sqrt(2*!pi*Sigma^2) * exp(-0.5*(ene-Emean)/Sigma^2)
		spec_type eq 'FileFunction'	: begin
						 	readcol,spec_file,logE,logdF
							if (min(logE) lt log_lemin) && (max(logE) gt log_lemax) then begin
							 spec=[]
							 for i=0,n_elements(logE)-2 do begin
							  w1 = where((lene ge logE(i)) and (lene lt logE(i+1)),w1_num)
							  if w1_num ne 0 then begin
							   pow1 = logdF(i) + (logdF(i+1)-logdF(i))/(logE(i+1)-logE(i))*(lene[w1]-logE(i))
							   spec=[spec,k*10^(pow1)]
							  endif
							 endfor
							endif else print, '### ERROR ### File"'+spec_file+'" does not cover total energy range.'
						 end
		else: print, '### ERROR ### wrong spec_type (ConstantValue/PowerLaw/PowerLaw2/ExpCutoff/BrokenPowerLaw/LogParabola/Gaussian/FileFunction)'
	endcase

	
	ebl_abs = fltarr(n_elements(ene))
        for m=0,n_elements(ene)-1 do   ebl_abs(m) =  exp(-tau_dominguez(ene(m)*1.d3,z_src))
	
	
	e_incr2 = (log_lemax-log_lemin)/(NEbins)
	lemin=[log_lemin:log_lemax-e_incr2:e_incr2]
	lemax=[log_lemin+e_incr2:log_lemax:e_incr2]
	flux=fltarr(n_elements(lemin))
	for i=0,n_elements(lemin)-1 do begin
		w=where(lene ge lemin(i) and lene lt lemax(i))
		flux(i)=meag(ebl_abs(w)*spec(w))*(10.^lemax(i)-10.^lemin(i))
	endfor

	plot_oo,ene,ebl_abs*spec,ytitle='dN / dE [ ph / cm!u2!n s TeV ]',xtitle='Photon Energy [TeV]'
	if spec_type eq 'FileFunction' then oplot,10^logE,k*10^logdF,psym=6,col=255 ; draw file_spec

endif else print, '### ERROR ### output file, (l,b) and (lemin,lemax) should be defined!'



file_delete, file, /allow_nonexistent

;;; >>>>>> Scrivi su file (if:diffuse/else:poit-like) <<<<<<

if keyword_set(DiffMap) then begin
	im = READFITS(DiffMap, im_header,EXTEN_NO=EXTEN_NO)
	

	map_size=size(im)
	sky=fltarr(map_size[1],map_size[2],n_elements(flux))
	pix_size1=abs(fxpar(im_header,'CDELT1'))
	pix_size2=abs(fxpar(im_header,'CDELT2'))
	for i=0,n_elements(flux)-1 do begin
		;sky[*,*,i]=im*flux[i] * pix_size1*pix_size2*(!pi/180)^2
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

