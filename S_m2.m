function [S_m] = S_m2(Alpha,X)
%   S_m calculate residual sum of squares
%   S_m = Sum(t=M+1:N)[X(N-t+1)-Alpha*X[N-t+2:N-t+M+1]]^2
%   where M is a order of approximation and N is quantity of points
M=length(Alpha);
N=length(X);

I = zeros(1,N-M);
Alpha1 = fliplr(Alpha);
N2 = N-M-1;


for i=1:M
    I = I + Alpha1(i)*X(i:N2+i);
end

% for t=1:N-M
% %     I(t-M)=(X(t)-X(t-M:t-1)*Alpha1)^2;
%     I(t)=X(t:t+M-1)*Alpha1;
% end

I = X(M+1:N) - I;

S_m=sum(I.^2);
