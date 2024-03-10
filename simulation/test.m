clear all;
close all;

format long;
tic;

%****************************  程序主体 *************
%%%%%%%%%%%%%%%%%  参数设定   %%%%%%%%%%%%%%
fc = 2000; %载波频率

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   信源    %%%%%%%%%%%%%%%%%%%
msg_source = [1 0 1 0 1 0 1 0 1 0 1 0 ];

figure(1);
subplot(3,1,1);
n = 0:length(msg_source)-1;
stem(n, msg_source, 'filled');
axis([0,12,-1,2]);  %画针状图
title('(a)信源时域波形');xlabel('n');ylabel('x(n)');
subplot(3,1,2);
Xk = fft(msg_source);
wn = 0:length(Xk)-1;
stem(2*wn./length(Xk), abs(Xk), 'filled');
axis([0,2,-1,10]);  %画针状图
title('(b)12点FFT得到的信源频域波形');xlabel('\omega/\pi');ylabel('X(k)');
subplot(3,1,3);
[Xk2,w]=freqz(msg_source, 1, 'whole'); % 'whole' 将区间设置为0-2pi, 否则是 0-pi。
plot(w/pi,abs(Xk2)); %创建 Y 中数据对 X 中对应值的二维线图
title('(c)freqz得到的信源频域波形');xlabel('\omega/\pi');ylabel('X(k)');


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












