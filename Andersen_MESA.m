function [Alpha,P] = Andersen_MESA( X )

%ANDRESEN_MESA uses Andersen algorithm to calculate parameters of the
% autoregresive model


   M=0;
   N=length(X);
   flX = fliplr(X);
   set(0,'RecursionLimit',N)
   FPE2=20000000000000;
%    c=0;
   FPE1=10000000000000;
%    A1=0;
%    A2=0; 
  [A1, c] = deal(0, 0);
  b1 = zeros(1,N);
  b2 = zeros(1,N);
  P=sum(X.^2)/N;
  
   while c<N/200 && c<1 && M<N/2    %with this cycle we didn't stop at local minimum of FPE
       if FPE2>FPE1
           A3=A1;
           c=0;
           FPE2=FPE1;
       else
           if M>N/4 %lowest order limit
               c=c+1;%check counter
               FPE2=FPE1;
           end
       end

       A2=A1;
       A1=0;

       M=M+1;

       if M==1
           b1(1)=X(1);
           b2(N-1) = X(N);
           b1(2:N-1) = X(2:N-1);
           b2(1:N-2) = X(2:N-1);
       else
%             A2_ = A2';
           b1(1:N-M) = b1_old(1:N-M)-A2(M-1)*b2_old(1:N-M);
           b2(1:N-M) = b2_old(2:N-M+1)-A2(M-1)*b1_old(2:N-M+1);

       end
       
       x=2*b1(1:N-M)*b2(1:N-M)'/(b1(1:N-M)*b1(1:N-M)'+b2(1:N-M)*b2(1:N-M)');
       b1_old=b1;
       b2_old=b2;
       if M==1
           A1=x;
           P=P*(1-A1^2);
       else
           A1=[A2-x*fliplr(A2), x];
           P=P*(1-A1(M)^2);
       end
              
       if M>N/4
           Sm=S_m2(A1,X);
           FPE1=FPE(N,M,Sm);
       end
   end
   %M=M-1;
   %A1=0;
   %for k=1:M+1
   %    p(k)=0;
   %   for t=1:N-k+1
   %          p(k)=p(k)+(X(t+k-1)-m)*(X(t)-m)/N;
   %    end
   %end
   %if M==1
   %    D=p(M+1)
   %    P=p(1)
   %    A1(M)=D/P
   %else
   %    D=p(M+1)-A3*fliplr(p(2:M))'
   %    P=p(1)-A3*p(2:M)'
   %    A1(M)=D/P
   %    A1(1:M-1)=A3-A1(M)*fliplr(A3)
   %end
   Alpha=A3';
   %Sm=S_m(Alpha,X);
   %FPE1=FPE(N,M,Sm)


end

