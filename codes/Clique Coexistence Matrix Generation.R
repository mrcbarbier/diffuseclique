### Clique diffuse coexistence matrix generation
   
MatrixGen=function(etai, beta.m, beta.sd)
{
S=length(etai)
mat=diag(S)
for(sp in 1:S)
	{ 
	Nmi=etai[-sp]
	rn=matrix(rnorm((S-2)^2),nrow=S-2, ncol=S-2)
	tmp=rbind(Nmi[-1],rn)
	AA=rbind(Nmi,t(tmp))
	QQ=t(qr.Q(qr(AA)))
	QQ=QQ*sign(QQ[1,1])
	x=QQ%*%rep(beta.m,S-1)
	x[1]=(1-etai[sp])/sqrt(sum(Nmi^2))
	x[2:(S-1)]=x[2:(S-1)]+beta.sd*rnorm(S-2)
	row=t(QQ)%*%x
	mat[sp,-sp]=row
	}
return(mat)
}


etai=c(0.2,0.6,0.8,0.9,0.1)
beta.m=0.5
beta.sd=0.4

m1=MatrixGen(etai, beta.m, beta.sd)
m1

