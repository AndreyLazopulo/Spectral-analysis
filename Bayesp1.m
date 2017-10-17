%find the postierior probability density (non-normalized) of parameters A,
%B when frequency is given from Bayesian analysis.
function y=Bayesp1(A,B,d,T,fs)
N=length(d);
t=0:1/fs:(N-1)/fs;
w=2*pi/T;
R=d*(cos(w*t))';
I=d*(sin(w*t))';
c=N/2+sin(N*w)/2/sin(w);
s=N/2-sin(N*w)/2/sin(w);
h1=R/sqrt(c);
h2=I/sqrt(s);
y=((sum(d.^2)-h1^2-h2^2+(A*sqrt(c)-h1).^2+(B*sqrt(s)-h2).^2)/N)^(-N/2);
