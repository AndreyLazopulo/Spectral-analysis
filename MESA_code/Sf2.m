function [f,Sf]=Sf2(Alpha,P,N,dt)

%SF2 calculates power spectrum Sf for given autoregression coeficients
%Alpha and varience of noise.

    if P==0
        fprintf('failed to calculate power spectrum\n')
    else
        i=sqrt(-1);
        T = (N-1)*dt;
        nf=2*N;
        
%         f_step=2/(4*N*dt);
%         T=2*dt:dt:50*60;
        M=length(Alpha);
        %alpha=[1:M];
        
        % calculate power spectrum
%         f=[0:f_step:1/(2*dt)]';
%         f=1./T';
        f=(1:nf)'./(T*4);
        term_j=exp(-i*2*pi*f*[1:M]*dt);
        term_k=term_j*Alpha;
        down=1.0-term_k;
        Sf=(P*dt)./(real(down).^2+imag(down).^2);
    end
end
