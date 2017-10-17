%logorithm of posterior of frequencies
function y=BayesLOGFRE(d,w,fs)
N=length(d);
t=-(N-1)/fs/2:1/fs:(N-1)/fs/2;
R=d*(cos(w*t))';
I=d*(sin(w*t))';
c=N/2+sin(N*w)/2/sin(w);
s=N/2-sin(N*w)/2/sin(w);
A=R/sqrt(c);
B=I/sqrt(s);
y=((2-N)/2)*log((1-(A^2+B^2)/sum(d.^2)));