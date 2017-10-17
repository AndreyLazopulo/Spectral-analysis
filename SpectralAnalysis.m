function [P,f] = SpectralAnalysis(samplingRate, Method)
%This function allows to calculate power spectra of data using one of four methods: MESA, FFT, Bayesian, Lomb-Scargle periodogram

%as output you write Sampling rate in seconds and for Method you write one
%of four methods: 'MESA', 'FFT', 'Bayesian' and 'LS' for Lomb-Scargle.

[filename,pathname]=uigetfile('.txt', 'Select simmulation files','MultiSelect','on');
    
    % find number of selected files
    if iscell(filename)
        NumberOfFiles=length(filename);
    elseif filename ~= 0
        NumberOfFiles = 1;
    else
        NumberOfFiles = 0;
    end
    
    
    for i=1:NumberOfFiles
    
        if NumberOfFiles == 1
            file1=fullfile(pathname,filename);
        else
            file1=fullfile(pathname,filename(i));
            file1=char(file1);
        end
        
         A=dlmread(file1);
         if size(A,1)>size(A,2)
            A=A';
         end
         
         switch Method
             case 'MESA'
                 samplingRate = samplingRate/60;
                 conf_MESA=@(N,p,v)exp((0.2285*log((1-p)/p)-1.325*p^0.14+0.3466*p^73.8+1.002)+(0.008666*(1-log((1-p)/p))+0.01156*log(p)+0.03183*p^1.874+0.1898)*log(N)+log(v))*samplingRate*3;
                 [P,f]=singleMESA(A,samplingRate,1,0,0);
                 figure;
                 plot(1./(f(1:end)*60),P(1:end));
                 f=f./60;
             case 'FFT'
                 L = length(A);
                 A=A-mean(A);
                 conf_FFT=@(N,p,v)exp(0.1187*log((1-p)./p)-0.5638*p.^0.3272+0.1465*p.^101.9+1.28).*N.^(0.005673*(1-log((1-p)./p))+0.004319*log(p)+0.02303*p.^1.151+0.1403)*v;
                 f_N=1/(2*samplingRate); % Determin Nyquist frequency
                 NFFT = 2^nextpow2(L); % Next power of 2 from length of y
                 Y1 = fft(A,NFFT);
                 P = abs(Y1).^2./NFFT;
                 f=0:2*f_N/length(P):2*f_N/length(P)*(length(P)-1);
                 figure;
                 plot(1./(f(1:length(P)/2)*3600),P(1:length(P)/2))
             case 'Bayesian'
                 samplingRate = 3600/samplingRate;
                 A=A-mean(A);
                 conf_Bayesian=@(N,p,V)exp(5.659*log((1-p)/p)+175.1*p^21.37+43.1*p^0.1864-45.25+(1.844*(1-log((1-p)/p))-0.9796*log(p)-6.146*p^0.3698+3.428)*log(N)+(0.06847*log((1-p)/p)+0.5182*p^0.2269-0.4491)*(log(N))^2+log(V))*180/samplingRate;
                 [~,P,f]=Bayesian2(A,samplingRate,4,4);
                 figure;
                 P = real(P);
                 plot(1./f,P)
                 f=f./3600;
             case 'LS'
                 t = (0:1:length(A)-1)*samplingRate;
                 conf_LS=@(N,p,V)exp((0.1261*log((1-p)/p)-0.4806*p^0.5805+0.1525*p^132.3+0.8118)+(0.009911*(1-log((1-p)/p))+0.002011*log(p)+0.05357*p^1.264+0.2614)*log(N)+(0.0003789*log((1-p)/p)-0.000634*p^1.382-0.001084*p^1.103-0.009852)*(log(N))^2+log(V));
                 [P,f]=fastlomb(A,t,0);
                 P=P*var(A);
                 figure;
                 plot(1./(f*3600),P)
             otherwise
                 P=NaN;
                 f=NaN;
                 error('No such method')
                 return
         end
    end
end

