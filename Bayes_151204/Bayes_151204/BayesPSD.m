function y=BayesPSD(d,w,fs,sigma,wmean,b)
N=length(d);
t=0:1/fs:(N-1)/fs;

r=length(wmean);
p=0;
for i=1:r
    C(i)=((d(i,:)*(cos(wmean(i)*t))')^2+(d(i,:)*(sin(wmean(i)*t))')^2)/N;
    p=p+(b(i)/2/pi/sigma)^(1/2)*exp(-b(i)*(wmean(i)-w)^2/2/sigma);
end
y=2*(sigma+sum(C))*p;

