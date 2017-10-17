data5=dlmread('clkAR_sleepbouts.txt');
data6=dlmread('clkJrk_sleepbouts.txt');
data4=dlmread('Hk1_fly1-100day5-15_LD_wakebouts.txt');
data3=dlmread('per0_sleepbouts.txt');

Nt=60;
fileID=fopen('hist_test_3.txt','a+');
for k=1:2
        fprintf(fileID,'per0 sleepbouts\n');
        tic
        [params,AIC,~,~]=PLwith3Expo_AL(data3,Nt,10);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        tic
        [params,AIC,outputstartvals,~,~]=PLwith3Expo_AL2(data3,20,0.05,20);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',outputstartvals);
        tic
        [params,AIC,~,~]=PLwith5Expo_AL(data3,Nt,10);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        tic
        [params,AIC,outputstartvals,~,~]=PLwith5Expo_AL2(data3,20,0.05,20);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',outputstartvals);
        fprintf(fileID,'Hk1_wakebouts\n');
        tic
        [params,AIC,~,~]=PLwith3Expo_AL(data4,Nt,10);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        tic
        [params,AIC,outputstartvals,~,~]=PLwith3Expo_AL2(data4,20,0.05,46);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',outputstartvals);
        tic
        [params,AIC,~,~]=PLwith5Expo_AL(data4,Nt,10);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        tic
        [params,AIC,outputstartvals,~,~]=PLwith5Expo_AL2(data4,Nt,0.05,46);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',outputstartvals);
        fprintf(fileID,'clkAR sleepbouts\n');
        tic
        [params,AIC,~,~]=PLwith3Expo_AL(data5,Nt,10);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        tic
        [params,AIC,outputstartvals,~,~]=PLwith3Expo_AL2(data5,Nt,0.05,20);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',outputstartvals);
        tic
        [params,AIC,~,~]=PLwith5Expo_AL(data5,Nt,10);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        tic
        [params,AIC,outputstartvals,~,~]=PLwith5Expo_AL2(data5,Nt,0.05,20);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',outputstartvals);
        fprintf(fileID,'clkJrk sleepbouts\n');
        tic
        [params,AIC,~,~]=PLwith3Expo_AL(data6,Nt,10);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        tic
        [params,AIC,outputstartvals,~,~]=PLwith3Expo_AL2(data6,Nt,0.05,20);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',outputstartvals);
        tic
        [params,AIC,~,~]=PLwith5Expo_AL(data6,Nt,10);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);  
        tic
        [params,AIC,outputstartvals,~,~]=PLwith5Expo_AL2(data6,Nt,0.05,20);
        toc
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',[params,AIC]);
        fprintf(fileID,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',outputstartvals);
end
fclose(fileID);
