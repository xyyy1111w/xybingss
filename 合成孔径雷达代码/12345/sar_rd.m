clc;
clear all;
close all;

%��������
Zc = 11e3;
Zd = 1.2e3;
Xmin =0;
Xmax =1.4e3;

H=10e3;

%��ʱ�����
c = 3e8;
fc = 9.6e9 ; %��Ƶ
lambda = c/fc;
v = 110 ;
% PRF = 233.9 ; %�����ظ�Ƶ��
% ds = 1/PRF ; %��ʱ��������
R0 = sqrt(Zc^2+H^2);
D = 1;
Lsar=lambda*R0/D ;  %�ϳɿ׾�����
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
Nslow=ceil((Xmax-Xmin+Lsar)/v/ds);%��������ceilΪȡ������
Nslow=2^nextpow2(Nslow);
% Nslow=4096;
s_min=(Xmin-Lsar/2)/v;
s_max=s_min+ds*(Nslow-1);
% s_max2=Lsar/2/v+DT;
s = linspace(s_min,s_max,Nslow);
fa=linspace(-0.5*PRF,0.5*PRF,Nslow);


%��ʱ�����

B = 300e6 ; %�������
tr = 2.5e-6 ; %����ʱ��
fs = 2*B ; %������
dt = 1/fs ;%��ʱ��������
K = B/tr ; %��Ƶ��
Rmin = sqrt((Zc-Zd).^2+H.^2);
Rmax = sqrt((Zc+Zd).^2+(Lsar/2).^2+H.^2); 
Nfast = ceil(2*(Rmax-Rmin)/c/dt+tr/dt); %��ʱ������� 
Nfast = 2^nextpow2(Nfast);
Nfast=Nfast/2;
% Nfast=11264;
t_min=2*Rmin/c;
t_max=t_min+dt*(Nfast-1);
t_max2=2*Rmax/c+tr;
t = linspace(t_min,t_max,Nfast);
fr=linspace(-0.5*fs,0.5*fs,Nfast);

Ntarget = 9; %Ŀ����Ŀ
target = [200 10e3 1
          200 11e3 1
          200 12e3 1
          700 10e3 1
          700 11e3 1
          700 12e3 1
          1200 10e3 1
          1200 11e3 1
          1200 12e3 1]; %Ŀ��λ��[x��z������ϵ��]

%���ɻز�
S0 = zeros(Nslow,Nfast);
for k=1:1:Ntarget
    sigma = target(k,3); %�õ�Ŀ��ķ���ϵ��
    Dslow = (s*v-target(k,1)); %��λ�����
    R = sqrt(Dslow.^2+target(k,2).^2+H.^2); %Ŀ��ʵ�ʾ���
    tau = 2*R/c; %ʱ���ӳ�
    Tfast = ones(Nslow,1)*t-tau'*ones(1,Nfast); %(t-tau),����չΪ����
    phase=pi*K*(Tfast).^2-(4*pi/lambda)*(R'*ones(1,Nfast)); %��λ
    S0=S0+sigma*exp(1i*phase).*(0<Tfast&Tfast<(tr)).*((abs(Dslow)<(Lsar/2))'*ones(1,Nfast));
   
end

% figure;
% imshow(real(S0));
clearvars phase Tfast

figure;  
imagesc(abs(S0));title('Amplitude of data') 
xlabel('������');ylabel('��λ��');
figure;
imagesc(angle(S0));title('Phase of data'); 
xlabel('������');ylabel('��λ��');

%������ѹ��
tu=t-2*Rmin/c;
Refr=exp(1i*pi*K*(tu).^2).*(0<tu&tu<(tr));  %�ο��źŲ���
% win=hanning((Nfast)).'; %����Ҫ�Ӵ�����ȥ��ע�ͣ���һ�м�ע��
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
xlabel('������'),ylabel('��λ��');
title('After Distance pulse compression');

%��λ����Ҷ�任
S_rd=ftx(Sr);
% S_rd=fftshift(fft(Sr,[],1),1);
% S_rd=zeros(Nslow,Nfast);
% for m=1:Nfast
%     S_rd(:,m)=fftshift(fft(Sr(:,m)));
% end

figure;
% colormap(gray)
imagesc(t*c/2,s*v,abs(S_rd));
xlabel('������'),ylabel('��λ��');
title('FFT of Azimuth');

% %���ξ���ѹ��
% tu=t-2*Rmin/c;
% Refr=exp(1i*pi*K*(tu).^2).*(0<tu&tu<(tr));  %�ο��źŲ���
% % win=hanning((Nfast)).'; %����Ҫ�Ӵ�����ȥ��ע�ͣ���һ�м�ע��
% win=1;
% Refr_conj=conj(fft((Refr.*win)));
% Sr=zeros(Nslow,Nfast);
% for M=1:1:Nslow
%     Sr(M,:)=ifft(fft(S_rd(M,:)).*Refr_conj);
% end
% % Sr=ifty(fty(S0).*(ones(Nslow,1)*conj(fty(Refr))));
% Gr=abs(Sr);


%-------------------------------------------------------------------------------
% %�����㶯У�� ������λ���
% tau_m=ones(Nfast,1)*t;
% f_eta_m=fa'*ones(1,Nslow);  %��λ��Ƶ�ʾ���
% f_tau_m=ones(Nfast,1)*fr;   %������Ƶ�ʾ���
% 
% % delta_R=lambda^2*fs/8/v^2*tau_m.*f_eta_m.^2;
% delta_R=lambda^2*R0/8/v^2*f_eta_m.^2;
% Grcmc=exp(1i*4*pi*f_tau_m.*delta_R/c);  %��λ�˷���
% temp=fftshift(fft(S_rd,[],2),2).*Grcmc;  %������λ���
% 
% clearvars Frcmc delta_R f_tau_m Grcmc
% clearvars f_tau_m
% 
% S_rcmc=ifftshift(ifft(temp,[],2),2);  %IFFT
%      
% clear temp
%------------------------------------------------------------------------------
%sinc��ֵ�����㶯У��
Kp = 1;
h = waitbar(0,'Sinc��ֵ');
P = 4;    %4��sinc��ֵ��
RMCmatrix = zeros(Nslow,Nfast);
for n=1:Nslow
    for m=P:Nfast
       delta_R = (1/8)*(lambda/v)^2*(R0+(m-Nfast/2)*c/2/fs)*((n-Nslow/2)*PRF/Nslow)^2;%���ȼ������Ǩ���� ���㷽�����ǰ�б��任��������������֪����
       RMC=2*delta_R*fs/c; %����ͽ���˼������뵥Ԫ
       delta_RMC = RMC-round(RMC); %����ͽ������С������
       for i = -P/2:P/2-1
           if m+RMC+i>Nfast %�ж��Ƿ񳬳��߽�
               RMCmatrix(n,m)=RMCmatrix(n,m)+S_rd(n,Nfast)*sinc(pi*(-i+RMC));
           else
               RMCmatrix(n,m)=RMCmatrix(n,m)+S_rd(n,m+round(RMC)+i)*sinc(pi*(-i+delta_RMC));
           end
       end
    end
    waitbar(n/Nslow);
end
close(h)

S_rmc = iftx((RMCmatrix));    %�����㶯У����ԭ��ʱ��
% S_rmc = ifft(fftshift(RMCmatrix,1),[],1);
Ga = abs(S_rmc);
%---------------------------------------------------------------------------------
figure;
% colormap(gray);
imagesc(t*c/2,s*v,abs(S_rmc));
xlabel('������');ylabel('��λ��');
title('Data After RCMC ');


%��λ��ѹ��
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
xlabel('������');ylabel('��λ��');
title('Final data');

Saup=upsample(Sa,4);

%-----------------------------------------------------------------------------------------------

plsr1=20*log10(abs(Sa(:,3460))/max(abs(Sa(:,3460)))); %Ŀ����ȸı�ʱ��Ҫ�޸�
figure;plot(plsr1);
xlabel('��λ��/m');ylabel('dB');
title('Pslr of azimuth direction');
% % plsr2=20*log10(abs(Sa(:,3460))/max(abs(Sa(:,3460)))); %Ŀ����ȸı�ʱ��Ҫ�޸�
% % figure;plot(plsr2);
% % xlabel('��λ��/m');ylabel('dB');
% % title('Pslr of azimuth direction');
% % plsr3=20*log10(abs(Sa(:,9034))/max(abs(Sa(:,9034)))); %Ŀ����ȸı�ʱ��Ҫ�޸�
% % figure;plot(plsr3);
% % xlabel('��λ��/m');ylabel('dB');
% % title('Pslr of azimuth direction');
plsr2=20*log10(abs(Sa(2866,:))/max(abs(Sa(2866,:)))); %Ŀ����ȸı�ʱ��Ҫ�޸�
figure;plot(plsr2);
xlabel('������');ylabel('dB');
title('Pslr of range direction');

















