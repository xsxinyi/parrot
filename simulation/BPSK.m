clear all;
close all;

format long;
tic;


%****************************  程序主体 *************
%%%%%%%%%%%%%%%%%  参数设定   %%%%%%%%%%%%%%
bit_rate = 1000;% 比特率
symbol_rate = 1000;%符号率
sps = 16;%每个符号的采样点数, sample per symbol
fc = 2000; %载波频率
fs = 16000; %采样频率

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   信源    %%%%%%%%%%%%%%%%%%%
%%%%% 随机信号 
msg_source = [ones(1,20) zeros(1,20) randi([0,1], 1, 99960)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   发射机  %%%%%%%%%%%%%%%%%%%%
bipolar_msg_source = 2*msg_source-1;

rollof_factor = 0.5;
rcos_fir = rcosdesign(rollof_factor, 6, sps);

%%插值
% for i=1:length(bipolar_msg_source)
%     up16_bipolar_msg_source(1+16*(i-1))=bipolar_msg_source(i);
%     up16_bipolar_msg_source(2+16*(i-1):16*i)=zeros(1,15);
% end
up16_bipolar_msg_source = upsample(bipolar_msg_source, 16);

rcos_msg_source = filter(rcos_fir, 1, up16_bipolar_msg_source);

% figure(1);
% plot(rcos_msg_source);
% title('时域波形');
% figure(2);
% plot(abs(fft(rcos_msg_source)));
% title('频域波形');

time = 1:length(rcos_msg_source);
rcos_msg_source_carrier = rcos_msg_source.*cos(2*pi*fc.*time/fs);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  信道  %%%%%%%%%%%%%%%%%%%%%%%%%
ebn0 = -6:8;
snr = ebn0 - 10*log10(0.5*16);

for i = 1:length(snr)
    rcos_msg_source_carrier_addnoise = awgn(rcos_msg_source_carrier, snr(i), "measured");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  接收机 %%%%%%%%%%%%%%%%%%%%%%%%%
    rcos_msg_source_addnoise = rcos_msg_source_carrier_addnoise.*cos(2*pi*fc.*time/fs);

    % figure(7)
    % plot(rcos_msg_source_addnoise);
    % title('时域波形');
    % figure(8);
    % plot(abs(fft(rcos_msg_source_addnoise)));
    % title('频域波形');

    fir_lp = fir1(128, 0.2);
    rcos_msg_source_lp = filter(fir_lp, 1, rcos_msg_source_addnoise);

    rollof_factor = 0.5;
    rcos_fir = rcosdesign(rollof_factor, 6, sps);

    rcos_msg_source_MF = filter(rcos_fir, 1, rcos_msg_source_lp);

    decision_site = 160;

    rcos_msg_source_MF_option = rcos_msg_source_MF(decision_site:sps:end);

    msg_source_MF_option_sign = sign(rcos_msg_source_MF_option);

    % figure(13);
    % plot(msg_source_MF_option_sign, '-*');
    % title('判决结果');
    % 
    % eyediagram(rcos_msg_source, sps);
    % title('发射端眼图');
    % eyediagram(rcos_msg_source_MF, sps);
    % title('接收端眼图');
    % 
    % scatterplot(rcos_msg_source(48+1:16:end-48));
    % title('BPSK星座图');

    %%%%%%%%%%%%%%%%%%%%%%%%% 信宿 %%%%%%%%%%%%%%%%%%%%%%%%%%%
    [err_number(i), bit_err_ratio(i)]=biterr(msg_source(1:length(rcos_msg_source_MF_option)), (msg_source_MF_option_sign+1)/2);
end
toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% 仿真结果 %%%%%%%%%%%%%%%%%%%%%%%%%
ber = berawgn(ebn0, 'psk', 2, 'nondiff');
semilogy(ebn0, bit_err_ratio, '-*', ebn0, ber, '-+');
xlabel('比特信噪比');
ylabel('误码率');
title('不同信噪比下误码率仿真曲线');
legend('实验曲线 ', '理论曲线')
grid on;












