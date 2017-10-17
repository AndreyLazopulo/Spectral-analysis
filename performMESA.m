function [P,f]=performMESA(samplingrate,binning,plotData,writeFile)

% PERFORMMESA function calculates power spectrum of the time series using
% Maximum Enthropy Spectral Analysis method (MESA) [1,2]. MESA uses autoregresive
% model as approximation for the time series. In this model datapoint at 
% time t is  presented as linear combination of M previous datapoint:
%
%                X(t)=a1*X(t-1)+a2*X(t-2)+...+aM*X(t-M),
%
% where M is the order of the approximation.
% This program uses Andersen algorithm [2] to calulate parameters a1-aM of
% the approximation. Order of approximation M is determined by minimizing
% the Final Prediction Error (FPE), calculated using the Akaike equation. 
% 
%SUBPROGRAMS
%   Andersen_MESA.m:  calculates autoregresive model parameters using
%                     Andersen algorhythm [2]
%   FPE.m:            calculates Final Prediction Error using Akaike
%                     equation [1]
%   Sf2.m:            calculates power spectrum from autoregresive model
%                     parameters [2]
%
%INPUTS
%   samplingrate:   sampling rate of the data in minutes
%   binning:        bin data to bigger time intervals (set 1 for no binning)
%   plotData:       if set to 1 the program will plot the data both in
%                   frequency domain and in time domain.
%   writeFile:      if set to 1 will save the power spectrum to txt file
%                   with freqquency in first column and power in second
%
%OUTPUTS
%   P:              MESA power spectrum
%   f:              respective frequencies (in 1/min)
%
%NOTES:
%   1. Data files should be single column txt files without header.
%   2. In order to simplify and accelerate calculations, the program
%   rescales data so that highest value in the time series is equal to 10. This
%   rescaling is accounted in power spectrum calculation, the obtained power
%   spectrum is scaled back in order to reflect the amplitude of the real data.
%
%REFERENCES
%   [1] Ulrych TJ, Bishop TN: Maximum entropy spectral analysis and autoregressive decomposition. Reviews of Geophysics 1975:183.
%   [2] Andersen N: On the calculation of filter coefficients for maximum entropy spectral analysis. Geophysics 1974:69�72.
%
%A. Lazopulo, Dec 2015

%% Select data files (one or more)
[filename,path1]=uigetfile('.txt', 'Select data files','MultiSelect','on');
file1=fullfile(path1,filename);
    % determin how many files you selected
    if iscell(filename)
        NumberOfFiles=length(filename);
    elseif filename ~= 0
        NumberOfFiles = 1;
        file1=cellstr(file1);
    else
        NumberOfFiles = 0;
    end
    
%% Calculate power spectrum for each file
    for i=1:NumberOfFiles
        Alpha3=0;
        X=0;
        Var2=0;
        
% 1) Open data file (file should be single column txt file)
       [~,name,~]=fileparts(char(file1(i)));
       graphname=name;
       disp(name)
       X=dlmread(char(file1(i)));
       
        if isempty(X)==1
            continue
        end
    
% 2) Bin data (if necessary)
        if binning~=1
            binner1=binning;
            timeseriestobin1 = X(1:(length(X)-mod(length(X),binner1)));
            binned1 = sum(reshape(timeseriestobin1,binner1,length(timeseriestobin1)/binner1));
            X = binned1';
        end
        
% 3) Data rescaling (see Notes part 2)
        X=X(5:end);
        if size(X,1)>size(X,2)
            X=X';
        end
        X=X-mean(X);
        scalingfactor=10/max(X);
        X=X*scalingfactor;
        
        
% 3) Perform MESA
        [Alpha3,Pm]=Andersen_MESA(X);
        [f,P]=Sf2(Alpha3,Pm,length(X),samplingrate*binning);
        disp('done')
        P=P./(scalingfactor^2); % scaling back the spectrum (see Notes part 2)
        
% 4) Plot power spectrum (if plotData == 1)
        if plotData == 1
            figure;
            subplot(2,1,1)
            plot(f(1:end),P(1:end))
            xlabel('frequency (1/min)')
            subplot(2,1,2)
            plot(1./(f(1:end)*60),P(1:end))
            xlim([0,40]);
            xlabel('period (hours)')
            title(graphname,'Interpreter', 'none')
            drawnow
        end
    
% 5) Save to file (if writeFile == 1)
        if writeFile == 1
            outputfile=[name,'_MESA_PS.txt'];
            File=[f,P];
            dlmwrite(outputfile,File,'\t') % write power spectrum in file
        end
    end