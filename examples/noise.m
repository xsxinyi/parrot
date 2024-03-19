clear all;
close all;

format long;


%**************************** 程序主体 *************
%%%%%%%%%%%%%%%%% 参数设定 %%%%%%%%%%%%%%
bit_rate = 1000;% 比特率
symbol_rate = 1000;%符号率
span = 6;
sps = 16;%每个符号的采样点数, sample per symbol
fc = 2000; %载波频率
fs = symbol_rate * sps; %采样频率
rollof_factor = 0.8;
rcos_fir = rcosdesign(rollof_factor, span, sps);
fir_M = 128;

msg_source = [1 0 0 1 1 0 1 1 1 0 1 1 1 0 1 0 1 1 0 0, randi([0,1], 1, 20000)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 发射机 %%%%%%%%%%%%%%%%%%%%
bipolar_msg_source = 2*msg_source-1;
bipolar_msg_source_I = bipolar_msg_source(1:2:end);
bipolar_msg_source_Q = bipolar_msg_source(2:2:end);
up16_bipolar_msg_source_I = upsample(bipolar_msg_source_I, sps);
up16_bipolar_msg_source_Q = upsample(bipolar_msg_source_Q, sps);
r_msg_source_I = filter(rcos_fir, 1, up16_bipolar_msg_source_I);
r_msg_source_Q = filter(rcos_fir, 1, up16_bipolar_msg_source_Q);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% 复数形式的加噪声 start %%%%%%%%%%%%%%%%%%
snr = 5;

time = 1:length(r_msg_source_I);
rcos_msg_source_carrier = (r_msg_source_I + 1j*r_msg_source_Q).*exp(1j*2*pi*fc.*time/fs);
%%%%% rcos_msg_source_carrier 的功率
spow1 = sum(rcos_msg_source_carrier*rcos_msg_source_carrier')/length(r_msg_source_I);
%%%线性高斯白噪声信道
%**加噪声
r_msg_source_carrier_addnoise = awgn(rcos_msg_source_carrier, snr ,'measured');
%%%%% rcos_msg_source_carrier 加噪声的功率
spow4 = sum(r_msg_source_carrier_addnoise*r_msg_source_carrier_addnoise')/length(r_msg_source_carrier_addnoise);
%%%% 噪声
noise_temp = r_msg_source_carrier_addnoise - rcos_msg_source_carrier;
%%%% 噪声功率
spow5 = sum(noise_temp*noise_temp')/length(noise_temp);
%%%% snr 
snr_temp = 10*log10(spow1/spow5);
%%%% 和 snr 相等
display(snr_temp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 复数形式的加噪声 end %%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%% 正交调制的加噪声 start %%%%%%%%%%%%%%%%

%****************************** 线性高斯白噪声信道
% *******方式一
%代码参考《simulation and software Radio for mobile communication》，书中是基带，这里是中频

EbN0 = 9;

%**************** QPSK 2bit, Es = 2Eb
k = 2;

tra_IFsignal = r_msg_source_I.*cos(2*pi*fc.*time/fs) - r_msg_source_Q.*sin(2*pi*fc.*time/fs);

%***
spow_S = sum(tra_IFsignal.^2)/length(tra_IFsignal);%中频信号功率
%***这个0.5是必须的, 因为假如单路的信号功率为1的话，这里正交调制的信号功率为2，A是根号2
%***S/N = Eb * k / N0 => kN/S = N0/Eb 
%***EbN0 = 10 * log10(Eb/N0)  => N0/Eb = 10.^(-EbN0/10)
%***N = N0/Eb = 10.^(-EbN0/10) * S/k 
%
attn_pow = sps * 0.5 * spow_S / k *10.^(-EbN0/10);  
attn = sqrt(attn_pow);
inoise = attn*randn(1,length(tra_IFsignal));
qnoise = attn*randn(1,length(tra_IFsignal));

IFsignal = tra_IFsignal+ inoise.*cos(2*pi*fc.*time/fs) - qnoise.*sin(2*pi*fc.*time/fs); %这句最关键

spow4 = sum(IFsignal*IFsignal')/length(IFsignal);
noise_temp = IFsignal - tra_IFsignal;
spow_N =  sum(noise_temp*noise_temp')/length(noise_temp);

snr_temp1 = 10*log10(spow_S/spow_N); %仿真代码中的信噪比
snr_temp2 = EbN0 +10*log10(k)-10*log10(0.5*fs/symbol_rate); %按照MATLAB的EbNO与SNR换算关系得到

display(snr_temp1);
display(snr_temp2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 正交调制的加噪声 end  %%%%%%%%%%%%%%%%%

