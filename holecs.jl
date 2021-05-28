using LinearAlgebra
L=41
t=-0.5
jk=1
jh=0.001*jk
T=0
h=zeros(L^2*2,L^2*2)
eta=10^-6
mu=-1.809
muf=0.000123
function theta(x,T=0)
	if T >0 
		fek=1/(1+exp(x/T))
	else 
		fek=x>0 ? 0 : 1
	end
return  fek
end


function get_order(U,ek,x,y,fz=1,T=0)
	
	if fz == 1
	result=sum([conj(U[x,i])*U[y,i]*(theta(ek[i],T)) for i in 1:length(ek)])
	
	else

	result=sum([U[x,i]*conj(U[y,i])*(1-theta(ek[i],T)) for i in 1:length(ek)])
 	end
 	return result

end

function post(x,L)
return (x[1]-1)*L+x[2]
end



function nearest(L)
	
n=[]

bh=[(i,j) for i in 1:L for j in 1:L]
p=[i for i in 1:L*L]

right=map(x->post(x,L),map(x->(x[2]>L ? (x[1],1) : x),map(x->x.+(0,1),bh)))

down=map(x->post(x,L),map(x->(x[1]>L ? (1,x[2]) : x),map(x->x.+(1,0),bh)))
a=map(tuple,p,right)
#println(a)
map(x->push!(n,x),a)
a=map(tuple,p,down)
map(x->push!(n,x),a)


return n
end 


function hole(n,L)
	centre=((L+1)/2,(L+1)/2)
centrepost=post(centre,L)
n=map(n) do x
	if sum(x.==centrepost) == 0
       x           
	end
end
deleteat!(n,n.==nothing)
return n 
end


nn=nearest(L)
nnhole=hole(nn,L)
c=[i for i in 1:2:(L)^2*2]
f=c.+1
cnext=[(c[i[1]],c[i[2]]) for i in nn]
fnext=[(f[i[1]],f[i[2]]) for i in nnhole]

for i in cnext
			h[i[1],i[2]]=t
			h[i[2],i[1]]=t
		end
#println(fnext)
#cdownnext=[(cup[i[1]],cup[i[2]]) for i in vcat(nn,cnn)]
#println(cupnext)
 


function hz(s,chi,muc,muf,L,n=1,T=0)
    MU=zeros(L*L*2,L*L*2)
    
    for i in c
		MU[i,i]=muc
		MU[i+1,i+1]=-muf
	end

	for kk in 1:n
    	global sget,chiget,nc,nf,errorc,errors,ss

	for i in 1:L*L
		h[c[i],f[i]]=-s[i]
		h[f[i],c[i]]=-s[i]
	end

	for i in 1:length(nnhole)
	
		h[fnext[i][1],fnext[i][2]]=-chi[i]
	
		h[fnext[i][2],fnext[i][1]]=-chi[i]

	end

	

	
	er,U=eigen(Hermitian(h-MU))
#println(er)
	nf=[get_order(U,er,i,i,1,T)  for i in f]
	nc=[get_order(U,er,i,i,1,T)  for i in c]
	sget=[jk*0.5*(get_order(U,er,i+1,i,1,T))  for i in c]
	ss=[jk*0.5*(get_order(U,er,i,i+1,1,T))  for i in c]
	chiget=[jh*0.5*(get_order(U,er,fnext[i][1],fnext[i][2],1,T))  for i in 1:length(nnhole)]
	errors=[abs((s[i]-sget[i])/(s[i]+im*eta)) for i in 1:length(s)] 
	errorc=[abs((chi[i]-chiget[i])/(chi[i]+im*eta)) for i in 1:length(chi)] 
    # if kk == n
    #     println(ss[errors.==maximum(errors)])
    #     println(s[errors.==maximum(errors)])
    #     println(sget[errors.==maximum(errors)])
    #      println(chi[errorc.==maximum(errorc)])
    #     println(chiget[errorc.==maximum(errorc)])
    #     println("----")
    # end
     s=sget
     chi=chiget
    end
return sget,chiget,nc,nf,errorc,errors,ss
end

 s,chi,nc,nf,errorc,errors,ss=hz(rand(L*L),rand(length(nnhole)),mu,muf,L,500)
fs=open("s.txt","w")
for i in s
    write(fs,string(i))
    write(fs,"\n")
end

close(fs)

fchi=open("chi.txt","w")
for i in chi
    write(fchi,string(i))
    write(fchi,"\n")
end

close(fchi)


fnc=open("nc.txt","w")
for i in nc
    write(fnc,string(i))
    write(fnc,"\n")
end

close(fnc)


fnf=open("nf.txt","w")
for i in nf
    write(fnf,string(i))
    write(fnf,"\n")
end

close(fnf)
ferrorc=open("errorc.txt","w")
for i in errorc
    write(ferrorc,string(i))
    write(ferrorc,"\n")
end

close(ferrorc)
ferrors=open("errors.txt","w")
for i in errors
    write(ferrors,string(i))
    write(ferrors,"\n")
end

close(ferrors)
fnn=open("nn.txt","w")
for i in nn
    write(fnn,string(i))
    write(fnn,"\n")
end

close(fnn)
