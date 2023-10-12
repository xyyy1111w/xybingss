clc;
clear all;
close all;

%成像区域
Zc = 11e3;
Zd = 1.2e3;
Xmin =0;
Xmax =1.4e3;

H=10e3;

%慢时间参数
c = 3e8;
fc = 9.6e9 ; %载频
lambda = c/fc;
v = 110 ;
% PRF = 233.9 ; %脉冲重复频率
% ds = 1/PRF ; %慢时间采样间隔
R0 = sqrt(Zc^2+H^2);
D = 1;
Lsar=lambda*R0/D ;  %合成孔径长度
Tsar=Lsar/v;

Ka=2*v^2/lambda/R0;
Ba=abs(Ka*Tsar);
PRF=Ba;
% PRF=233.9;
% Ka=PRF/Tsar;
PRT=1/PRF;
ds=PRT;

% DT=2;
% DX=DT*v;
Nslow=ceil((Xmax-Xmin+Lsar)/v/ds);%脉冲数，ceil为取整函数
Nslow=2^nextpow2(Nslow);
% Nslow=4096;
s_min=(Xmin-Lsar/2)/v;
s_max=s_min+ds*(Nslow-1);
% s_max2=Lsar/2/v+DT;
s = linspace(s_min,s_max,Nslow);
fa=linspace(-0.5*PRF,0.5*PRF,Nslow);


%快时间参数

B = 300e6 ; %发射带宽
tr = 2.5e-6 ; %发射时宽
fs = 2*B ; %采样率
dt = 1/fs ;%快时间采样间隔
K = B/tr ; %调频率
Rmin = sqrt((Zc-Zd).^2+H.^2);
Rmax = sqrt((Zc+Zd).^2+(Lsar/2).^2+H.^2); 
Nfast = ceil(2*(Rmax-Rmin)/c/dt+tr/dt); %快时间采样数 
Nfast = 2^nextpow2(Nfast);
Nfast=Nfast/2;
% Nfast=11264;
t_min=2*Rmin/c;
t_max=t_min+dt*(Nfast-1);
t_max2=2*Rmax/c+tr;
t = linspace(t_min,t_max,Nfast);
fr=linspace(-0.5*fs,0.5*fs,Nfast);

Ntarget = 9; %目标数目
target = [200 10e3 1
          200 11e3 1
          200 12e3 1
          700 10e3 1
          700 11e3 1
          700 12e3 1
          1200 10e3 1
          1200 11e3 1
          1200 12e3 1]; %目标位置[x，z，反射系数]

%生成回波
S0 = zeros(Nslow,Nfast);
for k=1:1:Ntarget
    sigma = target(k,3); %得到目标的反射系数
    Dslow = (s*v-target(k,1)); %方位向距离
    R = sqrt(Dslow.^2+target(k,2).^2+H.^2); %目标实际距离
    tau = 2*R/c; %时间延迟
    Tfast = ones(Nslow,1)*t-tau'*ones(1,Nfast); %(t-tau),并扩展为矩阵
    phase=pi*K*(Tfast).^2-(4*pi/lambda)*(R'*ones(1,Nfast)); %相位
    S0=S0+sigma*exp(1i*phase).*(0<Tfast&Tfast<(tr)).*((abs(Dslow)<(Lsar/2))'*ones(1,Nfast));
   
end

% figure;
% imshow(real(S0));
clearvars phase Tfast

figure;  
imagesc(abs(S0));title('Amplitude of data') 
xlabel('距离向');ylabel('方位向');
figure;
imagesc(angle(S0));title('Phase of data'); 
xlabel('距离向');ylabel('方位向');

%距离向压缩
tu=t-2*Rmin/c;
Refr=exp(1i*pi*K*(tu).^2).*(0<tu&tu<(tr));  %参考信号产生
% win=hanning((Nfast)).'; %若需要加窗，则去掉注释，下一行加注释
win=1;
Refr_conj=conj(fft((Refr.*win)));
Sr=zeros(Nslow,Nfast);
for M=1:1:Nslow
    Sr(M,:)=ifft(fft(S0(M,:)).*Refr_conj);
end
% Sr=ifty(fty(S0).*(ones(Nslow,1)*conj(fty(Refr))));
Gr=abs(Sr);

figure;
% colormap(gray)
imagesc((t)*c/2,s*v,Gr); 
xlabel('距离向'),ylabel('方位向');
title('After Distance pulse compression');

%方位向傅里叶变换
S_rd=ftx(Sr);
% S_rd=fftshift(fft(Sr,[],1),1);
% S_rd=zeros(Nslow,Nfast);
% for m=1:Nfast
%     S_rd(:,m)=fftshift(fft(Sr(:,m)));
% end

figure;
% colormap(gray)
imagesc(t*c/2,s*v,abs(S_rd));
xlabel('距离向'),ylabel('方位向');
title('FFT of Azimuth');

% %二次距离压缩
% tu=t-2*Rmin/c;
% Refr=exp(1i*pi*K*(tu).^2).*(0<tu&tu<(tr));  %参考信号产生
% % win=hanning((Nfast)).'; %若需要加窗，则去掉注释，下一行加注释
% win=1;
% Refr_conj=conj(fft((Refr.*win)));
% Sr=zeros(Nslow,Nfast);
% for M=1:1:Nslow
%     Sr(M,:)=ifft(fft(S_rd(M,:)).*Refr_conj);
% end
% % Sr=ifty(fty(S0).*(ones(Nslow,1)*conj(fty(Refr))));
% Gr=abs(Sr);


%-------------------------------------------------------------------------------
% %距离徙动校正 线性相位相乘
% tau_m=ones(Nfast,1)*t;
% f_eta_m=fa'*ones(1,Nslow);  %方位向频率矩阵
% f_tau_m=ones(Nfast,1)*fr;   %距离向频率矩阵
% 
% % delta_R=lambda^2*fs/8/v^2*tau_m.*f_eta_m.^2;
% delta_R=lambda^2*R0/8/v^2*f_eta_m.^2;
% Grcmc=exp(1i*4*pi*f_tau_m.*delta_R/c);  %相位乘法器
% temp=fftshift(fft(S_rd,[],2),2).*Grcmc;  %线性相位相乘
% 
% clearvars Frcmc delta_R f_tau_m Grcmc
% clearvars f_tau_m
% 
% S_rcmc=ifftshift(ifft(temp,[],2),2);  %IFFT
%      
% clear temp
%------------------------------------------------------------------------------
%sinc插值距离徙动校正
Kp = 1;
h = waitbar(0,'Sinc插值');
P = 4;    %4点sinc插值；
RMCmatrix = zeros(Nslow,Nfast);
for n=1:Nslow
    for m=P:Nfast
       delta_R = (1/8)*(lambda/v)^2*(R0+(m-Nfast/2)*c/2/fs)*((n-Nslow/2)*PRF/Nslow)^2;%首先计算距离迁移量 计算方法就是把斜距变换到距离多普勒域就知道了
       RMC=2*delta_R*fs/c; %距离徒动了几个距离单元
       delta_RMC = RMC-round(RMC); %距离徒动量的小数部分
       for i = -P/2:P/2-1
           if m+RMC+i>Nfast %判断是否超出边界
               RMCmatrix(n,m)=RMCmatrix(n,m)+S_rd(n,Nfast)*sinc(pi*(-i+RMC));
           else
               RMCmatrix(n,m)=RMCmatrix(n,m)+S_rd(n,m+round(RMC)+i)*sinc(pi*(-i+delta_RMC));
           end
       end
    end
    waitbar(n/Nslow);
end
close(h)

S_rmc = iftx((RMCmatrix));    %距离徙动校正后还原到时域
% S_rmc = ifft(fftshift(RMCmatrix,1),[],1);
Ga = abs(S_rmc);
%---------------------------------------------------------------------------------
figure;
% colormap(gray);
imagesc(t*c/2,s*v,abs(S_rmc));
xlabel('距离向');ylabel('方位向');
title('Data After RCMC ');


%方位向压缩
ta = s - (Xmin+791.7)/v;
Refa = exp(-1i*pi*Ka.*(ta).^2).*((-Tsar/2<ta&ta<Tsar/2));
win=hanning(Nslow).';
Sa=iftx(ftx(S_rmc).*(conj(ftx(Refa.*win)).'*ones(1,Nfast)));

% tau_m=ones(Nslow,1)*t;
% f_eta_m=fa'*ones(1,Nfast);
% Ka=2*v^2/lambda./(tau_m*c/2);
% Hf_ac=exp(-1i*pi*(f_eta_m).^2./Ka);
% Sac_f=S_rmc.*Hf_ac;
% Sa=iftx(Sac_f);

Gar = abs(Sa);

figure;
% colormap(gray);
imagesc(t*c/2,s*v,Gar);
xlabel('距离向');ylabel('方位向');
title('Final data');

Saup=upsample(Sa,4);

%-----------------------------------------------------------------------------------------------

plsr1=20*log10(abs(Sa(:,3460))/max(abs(Sa(:,3460)))); %目标深度改变时需要修改
figure;plot(plsr1);
xlabel('方位向/m');ylabel('dB');
title('Pslr of azimuth direction');
% % plsr2=20*log10(abs(Sa(:,3460))/max(abs(Sa(:,3460)))); %目标深度改变时需要修改
% % figure;plot(plsr2);
% % xlabel('方位向/m');ylabel('dB');
% % title('Pslr of azimuth direction');
% % plsr3=20*log10(abs(Sa(:,9034))/max(abs(Sa(:,9034)))); %目标深度改变时需要修改
% % figure;plot(plsr3);
% % xlabel('方位向/m');ylabel('dB');
% % title('Pslr of azimuth direction');
plsr2=20*log10(abs(Sa(2866,:))/max(abs(Sa(2866,:)))); %目标深度改变时需要修改
figure;plot(plsr2);
xlabel('距离向');ylabel('dB');
title('Pslr of range direction');

















