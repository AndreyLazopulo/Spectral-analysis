function performFFT(samplingrate,binning,plotData,writeFile)

[filename,path1]=uigetfile('.txt', 'Select data files','MultiSelect','on');

    % determin how many files you selected
    if iscell(filename)
        NumberOfFiles=length(filename);
    elseif filename ~= 0
        NumberOfFiles = 1;
    else
        NumberOfFiles = 0;
    end

    for i=1:NumberOfFiles
        Alpha3=0;
        X=0;
        Var2=0;
    
        if NumberOfFiles == 1
            filename
            file1=fullfile(path1,filename);
            [~,name,~]=fileparts(file1);
            graphname=filename;
        else
            filename(i)
            file1=fullfile(path1,filename(i));
            file1=char(file1);
            [~,name,~]=fileparts(file1);
            graphname=filename(i);
        end
        
        % open data file
       [~, X]=hdrload(file1);
       
        if isempty(X)==1
            continue
        end
    
        % bin data
        if binning~=1
            binner1=binning;
            timeseriestobin1 = X(1:(length(X)-mod(length(X),binner1)));
            binned1 = sum(reshape(timeseriestobin1,binner1,length(timeseriestobin1)/binner1));
            X = binned1';
        end
        
        L = length(X);
%           Fs=1/(samplingrate*binning); % Determin sampling rate
        NFFT = 2^nextpow2(L); % Next power of 2 from length of y
        Y1 = fft(X,NFFT);
        P = (abs(Y1(1:NFFT/2+1)).^2)./NFFT;
        f = 1/(2*samplingrate)*linspace(0,1,NFFT/2+1);
        if plotData==1
            figure;
            plot(1./(f*60),P)
        end
        
        if writeFile==1
            outputfile=[name,'_FFT_PS.txt'];
            File=[f',P];
            dlmwrite(outputfile,File,'\t') % write power spectrum in file
        end
    end
end