%logorithm of posterior of frequencies
function y=BayesLOGFRE2(d,w,fs)
N=length(d);
T=N/fs;
t=-(N-1)/fs/2:1/fs:(N-1)/fs/2;
R=d*(cos(w*t))';
I=d*(sin(w*t))';
c=N/2+sin(w*T)/2/sin(w*T/N);
s=N/2-sin(w*T)/2/sin(w*T/N);
A=R/sqrt(c);
B=I/sqrt(s);
y=((2-N)/2)*log((1-(A^2+B^2)/sum(d.^2)));