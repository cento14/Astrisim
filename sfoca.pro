function sfoca,l,b,lr,br,ins_file=ins_file,ene=ene,theta=theta

readcol,ins_file,lene_c,aeff_c,r68_c,r80_c,eres_c,bgrate_c,skip=1,/silent

ene_c=10.^lene_c

r68=interpol(r68_c/sqrt(2),ene_c,ene)
eres=interpol(eres_c,ene_c,ene)


lr = l + r68*randomu(seed, n_elements(l), /normal) /  cos(b/!radeg)
br = b + r68*randomu(seed, n_elements(b), /normal) 

ener = ene + eres*ene*randomu(seed, n_elements(ene), /normal) 

return,ener

end