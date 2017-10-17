% find the expectation value of A,B from posterior distribution when
% frequency is given from Bayesian analysis
function [A,B]=BayesPARA(d,T,fs)
N=length(d);
t=0:1/fs:(N-1)/fs;
w=2*pi/T;
R=d*(cos(w*t))';
I=d*(sin(w*t))';
c=N/2+sin(N*w)/2/sin(w);
s=N/2-sin(N*w)/2/sin(w);
A=R/c;
B=I/s;