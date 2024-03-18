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
fs = symbol_rate * sps; %采样频率
rollof_factor = 0.8;
rcos_fir = rcosdesign(rollof_factor, 6, sps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   信源    %%%%%%%%%%%%%%%%%%%
%%%%% 随机信号 
% msg_source = [ones(1,20) zeros(1,20) randi([0,1], 1, 99960)];
% msg_source = [ones(1,10) zeros(1,10) randi([0,1], 1, 10)];

msg_source = [1 0 0 1 1 0 1 1 1 0 1 1 1 0 1 0 1 1 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   发射机  %%%%%%%%%%%%%%%%%%%%
bipolar_msg_source = 2*msg_source-1;

figure(1);
stem(bipolar_msg_source, '.');axis([0,length(msg_source),-2,2]);
title('信源时域波形');
% figure(2);
% [H,w] = freqz(bipolar_msg_source, 1, 'whole');
% plot(w/pi,abs(H));
% title('信源freqz频域波形');

bipolar_msg_source_I = bipolar_msg_source(1:2:end);
bipolar_msg_source_Q = bipolar_msg_source(2:2:end);

% display(bipolar_msg_source_I);
% display(bipolar_msg_source_Q);

figure(3);
stem(bipolar_msg_source_I, '.');axis([0,length(bipolar_msg_source_I),-2,2]);
title('I路信源时域波形');
% figure(4);
% [H_I,w_I] = freqz(bipolar_msg_source_I, 1, 'whole');
% plot(w_I/pi,abs(H_I));
% title('I路信源freqz频域波形');

figure(5);
stem(bipolar_msg_source_Q, '.');axis([0,length(bipolar_msg_source_Q),-2,2]);
title('Q路信源时域波形');
% figure(6);
% [H_Q,w_Q] = freqz(bipolar_msg_source_Q, 1, 'whole');
% plot(w_Q/pi,abs(H_Q));
% title('Q路信源freqz频域波形');


%%%%%%%%%%%%% 上采样 %%%%%%%%%%%%%
up16_bipolar_msg_source_I = upsample(bipolar_msg_source_I, sps);
up16_bipolar_msg_source_Q = upsample(bipolar_msg_source_Q, sps);

% figure(7);
% stem(up16_bipolar_msg_source_I,'.');axis([0,length(up16_bipolar_msg_source_I),-2,2]);
% title('I路信源经插值后的时域波形');
% figure(8);
% [H_I2,w_I2] = freqz(up16_bipolar_msg_source_I, 1, 'whole');
% plot(w_I2/pi,abs(H_I2));
% title('I路信源经插值后freqz频域波形');
% figure(9);
% stem(up16_bipolar_msg_source_Q,'.');axis([0,length(up16_bipolar_msg_source_Q),-2,2]);
% title('Q路信源经插值后的时域波形');
% figure(10);
% [H_Q2,w_Q2] = freqz(up16_bipolar_msg_source_Q, 1, 'whole');
% plot(w_Q2/pi,abs(H_Q2));
% title('Q路信源经插值后freqz频域波形');

%%%%%%%%%%%%%%%% 成形 %%%%%%%%%%%%%%%
rcos_msg_source_I = filter(rcos_fir, 1, up16_bipolar_msg_source_I);
rcos_msg_source_Q = filter(rcos_fir, 1, up16_bipolar_msg_source_Q);

figure(11);
stem(rcos_msg_source_I,'.');
title('I路通过平方根升余弦成型滤波器后的时域波形');
% figure(12);
% [H_I4,w_I4] = freqz(rcos_msg_source_I, 1, 'whole');
% plot(w_I4/pi,abs(H_I4));
% title('I路通过平方根升余弦成型滤波器后的频域波形');
figure(13);
stem(rcos_msg_source_Q,'.');
title('Q路通过平方根升余弦成型滤波器后的时域波形');
% figure(14);
% [H_Q4,w_Q4] = freqz(rcos_msg_source_Q, 1, 'whole');
% plot(w_Q4/pi,abs(H_Q4));
% title('Q路通过平方根升余弦成型滤波器后的频域波形');


%%%%%载波发送
time = 1:length(rcos_msg_source_I);
tra_IFsignal = rcos_msg_source_I.*cos(2*pi*fc.*time/fs) - rcos_msg_source_I.*sin(2*pi*fc.*time/fs);

figure(15);
plot(tra_IFsignal);
title('上载波后的时域波形');
figure(16);
[H6,w6] = freqz(tra_IFsignal, 1, 'whole');
plot(w6/pi,abs(H6));
title('上载波后的频域波形');


%%%%%%相干解调
rcsv_msg_source_I = tra_IFsignal.*cos(2*pi*fc.*time/fs);
rcsv_msg_source_Q = tra_IFsignal.*(-sin(2*pi*fc.*time/fs));

figure(17)
plot(rcsv_msg_source_I);
% stem(rcos_msg_source_addnoise,'.');
title('相干解调后的I路时域波形');
figure(18);
[H7_I,w7_I] = freqz(rcsv_msg_source_I, 1, 'whole');
plot(w7_I/pi,abs(H7_I));
title('相干解调后的I路频域波形');

figure(19)
plot(rcsv_msg_source_Q);
% stem(rcos_msg_source_addnoise,'.');
title('相干解调后的Q路时域波形');
figure(20);
[H7_Q,w7_Q] = freqz(rcsv_msg_source_Q, 1, 'whole');
plot(w7_Q/pi,abs(H7_Q));
title('相干解调后的I路频域波形');

fir_lp = fir1(128, 0.2);
rcsv_msg_source_lp_I = filter(fir_lp, 1, rcsv_msg_source_I);
rcsv_msg_source_lp_Q = filter(fir_lp, 1, rcsv_msg_source_Q);

figure(21)
plot(rcsv_msg_source_lp_I);
% stem(rcos_msg_source_addnoise,'.');
title('FIR低通滤波后的I路时域波形');
figure(22);
[H8_I,w8_I] = freqz(rcsv_msg_source_lp_I, 1, 'whole');
plot(w8_I/pi,abs(H8_I));
title('FIR低通滤波后的I路频域波形');

figure(23)
plot(rcsv_msg_source_lp_Q);
% stem(rcos_msg_source_addnoise,'.');
title('FIR低通滤波后的Q路时域波形');
figure(24);
[H8_Q,w8_Q] = freqz(rcsv_msg_source_lp_Q, 1, 'whole');
plot(w8_Q/pi,abs(H8_Q));
title('FIR低通滤波后的I路频域波形');

rcsv_msg_source_MF_I = filter(rcos_fir, 1, rcsv_msg_source_lp_I);
rcsv_msg_source_MF_Q = filter(rcos_fir, 1, rcsv_msg_source_lp_Q);

figure(25)
% plot(rcsv_msg_source_MF_I);
stem(rcsv_msg_source_MF_I,'.');
title('均衡滤波后的I路时域波形');
figure(26);
[H9_I,w9_I] = freqz(rcsv_msg_source_MF_I, 1, 'whole');
plot(w9_I/pi,abs(H9_I));
title('均衡滤波后的I路频域波形');

figure(27)
% plot(rcsv_msg_source_MF_Q);
stem(rcsv_msg_source_MF_Q,'.');
title('均衡滤波后的Q路时域波形');
figure(28);
[H9_Q,w9_Q] = freqz(rcsv_msg_source_MF_Q, 1, 'whole');
plot(w9_Q/pi,abs(H9_Q));
title('均衡滤波后的I路频域波形');



