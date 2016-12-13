pro mksrc,file,l,b,lemin,lemax,flux,index=index,k=k,cutoff=cutoff


;  k =   flux @ 1 TeV    [ph / cm sec TeV]  

; flux = flusso integrato nel canale in [ph / cm sec]


if keyword_set(k) then begin

   if not keyword_set(index) then index=-2.
   if not keyword_set(cutoff) then cutoff=1e8


   lene=([findgen(3000)/1000.]-1.)
   ene=10.d0^lene
   spec=k*ene^index*exp(-ene/cutoff)

   lemin=findgen(300)*0.01-1. 
   lemax=lemin+0.01
   flux=fltarr(n_elements(lemin))

   for i=0,n_elements(lemin)-1 do begin
       
      w=where(lene ge lemin(i) and lene lt lemax(i))
      flux(i)=meag(spec(w))*(10.^lemax(i)-10.^lemin(i))

   endfor
   
   plot_oo,ene,spec,ytitle='dN / dE [ ph / cm!u2!n s TeV ]',xtitle='Photon Energy [TeV]'


endif



sky=fltarr(n_elements(l),n_elements(b),n_elements(flux))
sky(*)=flux

if n_elements(l) eq 1 then begin
   bin=1.
   endif else begin
         bin=l(1)-l(0)
         endelse

;>>>> Cancella eventuale file

;  if not keyword_set(ext) then begin
spawn,'\rm -rf '+file+' >& /dev/null' ; PR 2015-01-28
;  endif

;>>>>>> Scrivi su file

MKHDR,he1,sky,/image

;crval1=-l(n_elements(l)/2-2)
;crpix1=n_elements(l)/2-1
fxaddpar,he1,'CRVAL1',l(0),''
fxaddpar,he1,'CRPIX1',0.,''
fxaddpar,he1,'CDELT1',-bin,''
fxaddpar,he1,'CTYPE1','GLON-CAR',''
;fxaddpar,he1,'CTYPE1','GLON',''
fxaddpar,he1,'CUNIT1','deg     ',''

;fxaddpar,he1,'CRVAL2',b(n_elements(b)/2-2),''
;fxaddpar,he1,'CRPIX2',n_elements(b)/2-1,''
fxaddpar,he1,'CRVAL2',b(0),''
fxaddpar,he1,'CRPIX2',0,''
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
;fxaddpar,he1,'EQUINOX',2000,''
;fxaddpar,he1,'RADECSYS','FK5',''
fxaddpar,he1,'EXTNAME','FLUX'

writefits,file,sky,he1,/append

;;;;  Energy Bounds


FXBHMAKE,head,n_elements(lemin),'EBOUNDS'
;MKHDR,he2,ebou,/image
;fxaddpar,he2,'EXTNAME','EBOUNDS'

;fits_open,file,fcb,/append
;fits_write,fcb,[1.2,2],he2,extname='EBOUNDS'
;fits_close,fcb

fxbaddcol,1,head,float(1),'L_EMIN'
fxbaddcol,2,head,float(1),'L_EMAX'
fxbcreate,unit,file,head

  FOR I=0,n_elements(lemin)-1. do begin
  
  fxbwrite,unit,lemin(i),1,I+1
  fxbwrite,unit,lemax(i),2,I+1
  
  ENDFOR

fxbfinish,unit


end
