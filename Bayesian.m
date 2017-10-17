function [psd,f]=Bayesian(Data,Samplingrate,Order,Oversamplefactor)

fps=Samplingrate;
D=Data;
order=Order;
oaf=Oversamplefactor;
nof=length(Data);
t=[0:1:nof-1]/fps;

PP=polyfit(t,D,2);
Tr=polyval(PP,t);
d=D-Tr;%Detrend

w=2*pi*[1:oaf*nof/2]*fps/nof/oaf;
f=[1:oaf*nof/2]*fps/nof/oaf;
T=1./f;
prc=length(w);


res(1,:)=d;
for j=1:order
  for i=1:prc
    logy(j,i)=BayesLOGFRE2(res(j,:),w(i),fps);
  end

  logym=max(logy(j,:));
  logyn=logy(j,:)-logym;
  logT=log(sum(T*(exp(logyn))'))-log(sum(exp(logyn)));
  Tmean(j)=exp(logT);
  wmean(j)=2*pi/Tmean(j);
  wmean(j)=roundn(wmean(j),-3);
  logp(j,:)=logyn-log(sum(exp(logyn)));
  p(j,:)=exp(logp(j,:));
  [A(j),B(j)]=BayesPARA2(res(j,:),Tmean(j),fps);
  C(j)=sqrt(A(j)^2+B(j)^2);
  res(j+1,:)=res(j,:)-A(j)*cos(wmean(j)*t)-B(j)*sin(wmean(j)*t);
  b(j)=BayesPSDb(res(j,:),fps,wmean(j));
  v(j)=BayesVAR(res(j+1,:),A(j),B(j),fps);
end

for i=1:prc
   psd(i)=BayesPSD(res,w(i),fps,v(order),wmean(1:order),b);
end
figure;
plot(1./f,psd);






