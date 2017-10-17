function y=BayesLOGLIK(res,sigma)
N=length(res);
y=-N/2*log(2*pi)-N*log(sqrt(sigma))-sum(res.^2)/2/sigma;