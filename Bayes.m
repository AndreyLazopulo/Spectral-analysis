clear all;

fps=1;
order=1;

nof=10000;
t=[0:1:nof-1]/fps;
G=normrnd(1,1,size(t));

PP=polyfit(t,G,2);
Tr=polyval(PP,t);
d=G-Tr;%Detrend

oaf=1; 
w=2*pi*[1:oaf*nof/2]*fps/nof/oaf;
prc=length(w);
T=2*pi./w;


res(1,:)=d;
for j=1:order
  for i=1:prc
    logy(j,i)=BayesLOGFRE(res(j,:),w(i),fps);
  end

  logym=max(logy(j,:));
  logyn=logy(j,:)-logym;
  logT=log(sum(T*(exp(logyn))'))-log(sum(exp(logyn)));
  Tmean(j)=exp(logT);
  wmean(j)=2*pi/Tmean(j);
  wmean(j)=roundn(wmean(j),-3);
  logp(j,:)=logyn-log(sum(exp(logyn)));
  p(j,:)=exp(logp(j,:));
  [A(j),B(j)]=BayesPARA(res(j,:),Tmean(j),fps);
  C(j)=sqrt(A(j)^2+B(j)^2);
  res(j+1,:)=res(j,:)-A(j)*cos(wmean(j)*t)-B(j)*sin(wmean(j)*t);
  b(j)=BayesPSDb(res(j,:),fps,wmean(j));
  v(j)=BayesVAR(res(j+1,:),A(j),B(j),fps);
end

for i=1:prc
   psd(i)=BayesPSD(res,w(i),fps,v(order),wmean(1:order),b);
end

figure(100)
plot(T,psd);
title('DD Bayesian power spectral of grooming');xlabel('units:hours')





