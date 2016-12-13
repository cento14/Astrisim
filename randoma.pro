function randoma,nr,xd,yd,bin=bin


int=0
for i=0,n_elements(xd)-1 do  int=[int,int(-1)+yd(i)]
xint=[xd(0)-bin/2. , xd + bin/2. ]

p=randomu(ss,nr)
v=interpol(xint,int/total(yd),p)

return,v


end