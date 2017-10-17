function [f,Sf]=Sf2(Alpha,P,N,dt)

    %   This function plot power spectrum S_f for given autoregression coeficients
    %   Alpha and varience of noise.

    if P==0
        fprintf('failed to calculate power spectrum\n')
    else
        i=sqrt(-1);
        f_step=2/(4*N*dt);
%         T=2*dt:dt:50*60;
        M=length(Alpha);
        %alpha=[1:M];
        
        % calculate power spectrum
        f=[0:f_step:1/(2*dt)]';
%         f=1./T';
        term_j=exp(-i*2*pi*f*[1:M]*dt);
        term_k=term_j*Alpha;
        down=1.0-term_k;
        Sf=(P*dt)./(real(down).^2+imag(down).^2);
    end
end
