function [FPE]=FPE(N,M,S_m)
%FPE calculates final prediction error using Akaike equation:
%           FPE = S_m^2*(N+(M+1))/(N-(M+1))
	
	FPE = S_m^2*(N+(M+1))/(N-(M+1));