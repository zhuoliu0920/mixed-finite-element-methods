clear all
stage = 1;
n = 5;
[u,p,lambda] = pwcflargestage3_2(n);
[u1,p1,lambda1] =pwcflargestage3(n,0);
[e0,preal] = errorestimatep3(p1,n,0);
e1 = errorestimatep3(p,n,1);