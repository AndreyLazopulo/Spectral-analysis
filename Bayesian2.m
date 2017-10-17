function [wmean,psd,f]=Bayesian2(Data,Samplingrate,Order,Oversamplefactor)
%samplingrate = samples per hour;
%order = number of peaks looking for

fps=Samplingrate;
D=Data;
order=Order;
oaf=Oversamplefactor;
nof=length(Data);
Tmean=zeros(1,order);
wmean=zeros(1,order);
B1=zeros(1,order);
B2=zeros(1,order);
C0=zeros(1,order);
b=zeros(1,order);
v=zeros(1,order);

t=[0:1:nof-1]/fps;

% PP=polyfit(t,D,2);
% Tr=polyval(PP,t);
Tr=mean(D);
d=D-Tr;%Detrend
N=length(d);
% figure(1);
% hold on
% subplot(5,1,1)
% plot(t,d,'b')
% title('data')
w=2*pi*[1:oaf*nof/2]*fps/nof/oaf;
f=[1:oaf*nof/2]*fps/nof/oaf;
T=1./f;
prc=length(w);
logy=zeros(order,prc);

res(1,:)=d;

for j=1:order
    %     tic
    for i=1:prc
        logy(j,i)=BayesLOGFRE2(res(j,:),w(i),fps);
    end
%     toc
%     figure;
%     plot(T,logy(j,:))
    if j>1
        for k=1:j-1
            [~,ind]=min(abs(w-wmean(k)));
            [~,ind2]=findpeaks(logy(1,:));
            if isempty(ind2(find(ind2-ind<0,1,'last')))==0
                logy(j,ceil((ind2(find(ind2-ind<0,1,'last'))+ind2(find(ind2-ind==0,1)))/2):floor((ind2(find(ind2-ind>0,1,'first'))-2+ind2(find(ind2-ind==0,1)))/2))=0;
            else
                logy(j,1:floor((ind2(find(ind2-ind>0,1,'first'))-2+ind2(find(ind2-ind==0,1)))/2))=0;
            end
        end
    end
    logym=max(logy(j,:));
    logyn=logy(j,:)-logym;
%     figure;
%     plot(T,logy(j,:));
    logT=log(sum(T*(exp(logyn))'))-log(sum(exp(logyn)));
    Tmean(j)=exp(logT);
    wmean(j)=2*pi/Tmean(j);
%     [~,ind]=min(abs(w-wmean(j)))
    wmean(j)=roundn(wmean(j),-3);
    [B1(j),B2(j)]=BayesPARA2(res(j,:),Tmean(j),fps);
    %   C(j)=sqrt(B1(j)^2+B2(j)^2);
    %   subplot(5,1,j);
    %   figure(j);
    %   hold on
    %   plot(t,B1(j)*cos(wmean(j)*t)+B2(j)*sin(wmean(j)*t),'r')
    res(j+1,:)=res(j,:)-B1(j)*cos(wmean(j)*t)-B2(j)*sin(wmean(j)*t);
    %   subplot(5,1,j+1)
    %   figure(j+1);
    %   plot(t,res(j+1,:))
    %   title('after substraction')
    v(j)=BayesVAR(res(j+1,:),B1(j),B2(j),fps);
    %   w(ind-10:ind+10)=0;
end
v_end=v(end);
p=0;
for i=1:order
    C0(i)=((d*(cos(wmean(i)*t))')^2+(d*(sin(wmean(i)*t))')^2)/N;
    %     b(i)=-B(i)^2*((sum(d.*t.*cos(wmean(i)*t)))^2-sum(d.*sin(wmean(i)*t))*sum(d.*(t.^2).*sin(wmean(i)*t)))-A(i)^2*((sum(d.*t.*sin(wmean(i)*t)))^2-sum(d.*cos(wmean(i)*t))*sum(d.*(t.^2).*cos(wmean(i)*t)));
    A=sum(d.*t.*cos(wmean(i)*t));
    B=sum(d.*t.*sin(wmean(i)*t));
    C=sum(d.*sin(wmean(i)*t))*sum(d.*(t.^2).*sin(wmean(i)*t));
    D=sum(d.*cos(wmean(i)*t))*sum(d.*(t.^2).*cos(wmean(i)*t));
    %     A=(d.*t)*(sin(wmean(i)*t))';
    %     B=(d.*t)*(cos(wmean(i)*t))';
    %     C=(d*(cos(wmean(i)*t))')*((d.*t.*t)*(cos(wmean(i)*t))');
    %     D=(d*(sin(wmean(i)*t))')*((d.*t.*t)*(sin(wmean(i)*t))');
    b(i)=-(A^2+B^2-C-D)*2/(N);
    %     term(i,:)=(b(i)/(2*pi*v))^0.5.*exp(-(b(i))*(wmean(i)-w).^2/(2*v));
    if b(i)>0
        if i>1 && abs(v(i)-v(i-1))<v(i-1)*0.008
            warning(['the peak at ',num2str(2*pi/wmean(i)),' is likely produced by noise']);
        end
        p=p+exp(log((b(i)/2/pi/v_end))*(1/2)+(-b(i)*(wmean(i)-w).^2/2/v_end));
    end
end
% disp(2*pi./wmean)
% disp(C0)
% disp(b)
% disp(v_end)
psd=2*(v_end+sum(C0))*p;
% figure;
% plot(2*pi./w,psd);
end



