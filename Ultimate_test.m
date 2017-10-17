T=28;
t=0:1:15000;
% samplingRate=1;
conf_LS=@(N,p,V)exp(0.6075*exp(-45.33*p)+1.008+(0.01344*log(p)+0.2773)*log(N)+(0.0002174*p^-0.3791-0.009472)*(log(N))^2+log(V));
conf_MESA=@(N,p,v)exp((1.097*exp(-135.2*p)-6.038*p+1.132)+(0.02172*log(p)+0.2036)*log(N)+log(v))*1*3;
conf_FFT=@(N,p,v)exp(0.497*exp(-115.3*p)-3.159*p+1.584)*N^(0.01088*log(p)+0.1496)*v;
conf_Bayesian=@(N,p,V)exp(6.418*exp(-89.08*p)-15.17*p-3.235+(-0.05523*(log(p))^2-2.328*p+1.36)*log(N)+(0.00327*(log(p))^2+0.1385*p-0.02145)*(log(N))^2+log(V))*180/60;
for i=1:27
    %MESA
    X=square(2*pi*t/((T-i)*60))+normrnd(0,4,size(t));
    tic
    [P,f]=singleMESA(X,1,1,0,0);
    t_MESA(i)=toc;
    [Pmax,ind]=max(P);
    acc_MESA(i)=(1/(60*f(ind))-(T-i))/(T-i);
    sens_MESA(i)=Pmax/conf_MESA(10000,0.05,var(X));
    %Lomb-Scargle
    tic
    [P,f]=fastlomb(X,t.*60,0);
    t_LS(i)=toc;
    P=P*var(X);
    [Pmax,ind]=max(P);
    acc_LS(i)=(1/(3600*f(ind))-(T-i))/(T-i);
    sens_LS(i)=Pmax/conf_LS(10000,0.05,var(X));
    %DFT
    L = length(X);
    X=X-mean(X);
    tic
    f_N=1/2; % Determin Nyquist frequency
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y1 = fft(X,NFFT);
    P = abs(Y1).^2./NFFT;
    f=0:2*f_N/length(P):2*f_N/length(P)*(length(P)-1);
    t_FFT(i)=toc;
    [Pmax,ind]=max(P);
    acc_FFT(i)=(1/(60*f(ind))-(T-i))/(T-i);
    sens_FFT(i)=Pmax/conf_FFT(10000,0.05,var(X));
    %Bayesian
%     X=X-mean(X);
    tic
    [~,P,f]=Bayesian2(X,60,4,4);
    t_Bayes(i)=toc;
    [Pmax,ind]=max(P);
    acc_Bayes(i)=(2*pi/(f(ind))-(T-i))/(T-i);
    sens_Bayes(i)=Pmax/conf_Bayesian(10000,0.05,var(X));
end
    