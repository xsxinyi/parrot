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
rollof_factor = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   信源    %%%%%%%%%%%%%%%%%%%
%%%%% 随机信号 
msg_source = [ones(1,10) zeros(1,10) randi([0,1], 1, 10)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%   发射机  %%%%%%%%%%%%%%%%%%%%
bipolar_msg_source = 2*msg_source-1;

figure(1);
stem(bipolar_msg_source, '.');axis([0,29,-2,2]);
title('信源时域波形');
% figure(2);
% n=0:29;
% stem(n/30*2, abs(fft(bipolar_msg_source)), '.');
% title('信源fft频域波形');
% figure(3);
% [H,w] = freqz(bipolar_msg_source, 1, 'whole');
% plot(w/pi,abs(H));
% title('信源freqz频域波形');

rcos_fir = rcosdesign(rollof_factor, 6, sps);

%%插值
up16_bipolar_msg_source = upsample(bipolar_msg_source, sps);

figure(4);
stem(up16_bipolar_msg_source,'.');axis([0,30*sps,-2,2]);
title('信源经插值后的时域波形');
% figure(5);
% [H2,w2] = freqz(up16_bipolar_msg_source, 1, 'whole');
% plot(w2/pi,abs(H2));
% title('信源经插值后freqz频域波形');

% figure(6);
% stem(rcos_fir,'.');axis([0,sps*6,-2,2]);
% title('平方根升余弦滤波器的时域波形');
% figure(7);
% [H3,w3] = freqz(rcos_fir, 1, 'whole');
% plot(w3/pi,20*log10(abs(H3)));
% title('平方根升余弦滤波器的频域波形');

rcos_msg_source = filter(rcos_fir, 1, up16_bipolar_msg_source);

figure(8);
% plot(rcos_msg_source);
stem(rcos_msg_source,'.');
title('通过平方根升余弦滤波器后的时域波形');
% figure(9);
% [H4,w4] = freqz(rcos_msg_source, 1, 'whole');
% plot(w4/pi,abs(H4));
% % plot(abs(fft(rcos_msg_source)));
% % stem(abs(fft(rcos_msg_source)),'.');
% %%%%%%%%% 可以看出截了最初的信源频谱的0~pi部分。
% title('通过平方根升余弦滤波器后的频域波形');

time = 1:length(rcos_msg_source);
rcos_msg_source_carrier = rcos_msg_source.*cos(2*pi*fc.*time/fs);



% figure(10);
% % stem(time, rcos_msg_source.*cos(2*pi*fc.*time/fs),'.');axis([-1,128,-0.4,0.4]);
% stem(rcos_msg_source_carrier,'.');
% title('上载波后的时域波形');
% figure(11);
% [H5,w5] = freqz(rcos_msg_source_carrier, 1, 'whole');
% plot(w5/pi,abs(H5));
% title('上载波后的频域波形');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  信道  %%%%%%%%%%%%%%%%%%%%%%%%%
ebn0 = -6:8;
snr = ebn0 - 10*log10(0.5*16);

% for i = 1:length(snr)
    %%%% awgn(x,snr,'measured')，首先计算输入x信号的功率，按照snr添加相应功率的高斯白噪声。
    % rcos_msg_source_carrier_addnoise = awgn(rcos_msg_source_ca
    % rrier, snr(i), "measured");
    rcos_msg_source_carrier_addnoise = awgn(rcos_msg_source_carrier, -6, "measured");

    figure(12);
    % stem(rcos_msg_source_carrier_addnoise,'.');
    plot(rcos_msg_source_carrier_addnoise);
    title('上载波加噪声后的时域波形');
    % figure(13);
    % [H6,w6] = freqz(rcos_msg_source_carrier_addnoise, 1, 'whole');
    % plot(w6/pi,abs(H6));
    % title('上载波加噪声后的频域波形');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  接收机 %%%%%%%%%%%%%%%%%%%%%%%%%
    rcos_msg_source_addnoise = rcos_msg_source_carrier_addnoise.*cos(2*pi*fc.*time/fs);

    figure(14)
    plot(rcos_msg_source_addnoise);
    % stem(rcos_msg_source_addnoise,'.');
    title('相干解调后的接收时域波形');
    % figure(15);
    % [H7,w7] = freqz(rcos_msg_source_addnoise, 1, 'whole');
    % plot(w7/pi,abs(H7));
    % title('相干解调后的接收频域波形');


    fir_lp = fir1(128, 0.2);
    rcos_msg_source_lp = filter(fir_lp, 1, rcos_msg_source_addnoise);

    figure(16)
    % plot(rcos_msg_source_addnoise);
    stem(rcos_msg_source_lp,'.');
    title('经过低通FIR滤波器后的接收时域波形');
    % figure(17);
    % [H8,w8] = freqz(rcos_msg_source_lp, 1, 'whole');
    % plot(w8/pi,abs(H8));
    % title('经过低通FIR滤波器后的接收频域波形');

    rcos_fir = rcosdesign(rollof_factor, 6, sps);
    rcos_msg_source_MF = filter(rcos_fir, 1, rcos_msg_source_lp);

    figure(18)
    % plot(rcos_msg_source_addnoise);
    stem(rcos_msg_source_MF,'.');
    title('经过平方根余弦滤波器后的接收时域波形');
    % figure(19);
    % [H9,w9] = freqz(rcos_msg_source_MF, 1, 'whole');
    % plot(w9/pi,abs(H9));
    % title('经过平方根余弦滤波器后的接收频域波形');

    decision_site = 160;
    rcos_msg_source_MF_option = rcos_msg_source_MF(decision_site:sps:end);

    figure(20)
    % plot(rcos_msg_source_addnoise);
    stem(rcos_msg_source_MF_option,'.');
    title('经过延时处理并抽取后的接收时域波形');
    % figure(21);
    % [H10,w10] = freqz(rcos_msg_source_MF_option, 1, 'whole');
    % plot(w10/pi,abs(H10));
    % title('经过延时处理并抽取后的接收频域波形');

    msg_source_MF_option_sign = sign(rcos_msg_source_MF_option);

    figure(22);
    plot(msg_source_MF_option_sign, '-*');
    title('判决结果');
    % 
    % eyediagram(rcos_msg_source, sps);
    % title('发射端眼图');
    % eyediagram(rcos_msg_source_MF, sps);
    % title('接收端眼图');
    % 
    % scatterplot(rcos_msg_source(48+1:16:end-48));
    % title('BPSK星座图');

    %%%%%%%%%%%%%%%%%%%%%%%%% 信宿 %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[err_number(i), bit_err_ratio(i)]=biterr(msg_source(1:length(rcos_msg_source_MF_option)), (msg_source_MF_option_sign+1)/2);
    bit_err_ratio = biterr(msg_source(1:length(rcos_msg_source_MF_option)), (msg_source_MF_option_sign+1)/2);
% end
toc;

display(bit_err_ratio);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% 仿真结果 %%%%%%%%%%%%%%%%%%%%%%%%%
% ber = berawgn(ebn0, 'psk', 2, 'nondiff');
% semilogy(ebn0, bit_err_ratio, '-*', ebn0, ber, '-+');
% xlabel('比特信噪比');
% ylabel('误码率');
% title('不同信噪比下误码率仿真曲线');
% legend('实验曲线 ', '理论曲线')
% grid on;












