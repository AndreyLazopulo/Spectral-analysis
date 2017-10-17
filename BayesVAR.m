function y=BayesVAR(d,A,B,fs)
N=length(d);
t=0:1/fs:(N-1)/fs;
m=length(A)+length(B);
y=(sum(d.^2)-sum(A.^2)-sum(B.^2))/(N-m-2);