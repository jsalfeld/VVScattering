BEGIN {
sum1=0
sum2=0
sum3=0
ntot=0
name="dummy"
}
{
sum1=sum1+$1/($2*$2)
sum3=sum3+$1*$1/($2*$2)
sum2=sum2+1./($2*$2)
ntot=ntot+1
name=$3
}
END{
res=sum1/sum2
del=1/sqrt(sum2)
chi=sum3-sum1**2/sum2
if(chi<0) chi=0;
printf("%100s: %20.15f +/- %20.15f / %5.2f\n",name,res,del,chi/ntot);
}
