clear global;
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
span = 6;
rcos_fir = rcosdesign(rollof_factor, span, sps);
fir_M = 128;
m = 2; % QPSK,一个符号2bit.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   信源    %%%%%%%%%%%%%%%%%%%
%%%%% 随机信号 
msg_source = [ones(1,20) zeros(1,20) randi([0,1], 1, 999960)];
% msg_source = [ones(1,10) zeros(1,10) randi([0,1], 1, 10)];

% msg_source = [1 0 0 1 1 0 1 1 1 0 1 1 1 0 1 0 1 1 0 0, randi([0,1], 1, 20)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   发射机  %%%%%%%%%%%%%%%%%%%%
bipolar_msg_source = 2*msg_source-1;

% figure(1);
% stem(bipolar_msg_source, '.');axis([0,length(msg_source),-2,2]);
% title('信源时域波形');
% figure(2);
% [H,w] = freqz(bipolar_msg_source, 1, 'whole');
% plot(w/pi,abs(H));
% title('信源freqz频域波形');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   I路和Q路  %%%%%%%%%%%%%%%%
bipolar_msg_source_I = bipolar_msg_source(1:2:end);
bipolar_msg_source_Q = bipolar_msg_source(2:2:end);

% display(bipolar_msg_source_I);
% display(bipolar_msg_source_Q);

% figure(3);
% stem(bipolar_msg_source_I, '.');axis([0,length(bipolar_msg_source_I),-2,2]);
% title('I路信源时域波形');
% figure(4);
% [H_I,w_I] = freqz(bipolar_msg_source_I, 1, 'whole');
% plot(w_I/pi,abs(H_I));
% title('I路信源freqz频域波形');

% figure(5);
% stem(bipolar_msg_source_Q, '.');axis([0,length(bipolar_msg_source_Q),-2,2]);
% title('Q路信源时域波形');
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

% figure(11);
% stem(rcos_msg_source_I,'.');
% title('I路通过平方根升余弦成型滤波器后的时域波形');
% figure(12);
% [H_I4,w_I4] = freqz(rcos_msg_source_I, 1, 'whole');
% plot(w_I4/pi,abs(H_I4));
% title('I路通过平方根升余弦成型滤波器后的频域波形');
% figure(13);
% stem(rcos_msg_source_Q,'.');
% title('Q路通过平方根升余弦成型滤波器后的时域波形');
% figure(14);
% [H_Q4,w_Q4] = freqz(rcos_msg_source_Q, 1, 'whole');
% plot(w_Q4/pi,abs(H_Q4));
% title('Q路通过平方根升余弦成型滤波器后的频域波形');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  信道  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%载波发送
time = 1:length(rcos_msg_source_I);
tra_IFsignal = rcos_msg_source_I.*cos(2*pi*fc.*time/fs) - rcos_msg_source_Q.*sin(2*pi*fc.*time/fs);

% figure(15);
% plot(tra_IFsignal);
% title('上载波后的时域波形');
% figure(16);
% [H6,w6] = freqz(tra_IFsignal, 1, 'whole');
% plot(w6/pi,abs(H6));
% title('上载波后的频域波形');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 加噪声 %%%%%%%%%%%%%%%%%%%%%%%%
ebn0 = 3:8;
spow_S = sum(tra_IFsignal.^2)/length(tra_IFsignal);%中频信号功率

err_number = zeros(1, length(ebn0));
bit_err_ratio = zeros(1, length(ebn0));

for j = 1:length(ebn0)
    % 实数的白噪声功率谱密度是 n0/2, 所以乘以 0.5,
    attn_pow = sps * 0.5 * spow_S / m *10.^(-ebn0(j)/10);  
    attn = sqrt(attn_pow);
    inoise = attn*randn(1,length(tra_IFsignal));
    qnoise = attn*randn(1,length(tra_IFsignal));

    tra_IFsignal_add_noise = tra_IFsignal+ inoise.*cos(2*pi*fc.*time/fs) - qnoise.*sin(2*pi*fc.*time/fs); %这句最关键


    %%%%%%相干解调
    rcsv_msg_source_I = tra_IFsignal_add_noise.*cos(2*pi*fc.*time/fs);
    rcsv_msg_source_Q = tra_IFsignal_add_noise.*(-sin(2*pi*fc.*time/fs));

    % figure(17)
    % plot(rcsv_msg_source_I);
    % % stem(rcos_msg_source_addnoise,'.');
    % title('相干解调后的I路时域波形');
    % figure(18);
    % [H7_I,w7_I] = freqz(rcsv_msg_source_I, 1, 'whole');
    % plot(w7_I/pi,abs(H7_I));
    % title('相干解调后的I路频域波形');

    % figure(19)
    % plot(rcsv_msg_source_Q);
    % % stem(rcos_msg_source_addnoise,'.');
    % title('相干解调后的Q路时域波形');
    % figure(20);
    % [H7_Q,w7_Q] = freqz(rcsv_msg_source_Q, 1, 'whole');
    % plot(w7_Q/pi,abs(H7_Q));
    % title('相干解调后的I路频域波形');

    fir_lp = fir1(fir_M, 0.2);
    rcsv_msg_source_lp_I = filter(fir_lp, 1, rcsv_msg_source_I);
    rcsv_msg_source_lp_Q = filter(fir_lp, 1, rcsv_msg_source_Q);

    % figure(21)
    % % plot(rcsv_msg_source_lp_I);
    % stem(rcsv_msg_source_lp_I,'.');
    % title('FIR低通滤波后的I路时域波形');
    % figure(22);
    % [H8_I,w8_I] = freqz(rcsv_msg_source_lp_I, 1, 'whole');
    % plot(w8_I/pi,abs(H8_I));
    % title('FIR低通滤波后的I路频域波形');

    % figure(23)
    % % plot(rcsv_msg_source_lp_Q);
    % stem(rcsv_msg_source_lp_Q,'.');
    % title('FIR低通滤波后的Q路时域波形');
    % figure(24);
    % [H8_Q,w8_Q] = freqz(rcsv_msg_source_lp_Q, 1, 'whole');
    % plot(w8_Q/pi,abs(H8_Q));
    % title('FIR低通滤波后的I路频域波形');

    rcsv_msg_source_MF_I = filter(rcos_fir, 1, rcsv_msg_source_lp_I);
    rcsv_msg_source_MF_Q = filter(rcos_fir, 1, rcsv_msg_source_lp_Q);

    % figure(25)
    % % plot(rcsv_msg_source_MF_I);
    % stem(rcsv_msg_source_MF_I,'.');
    % title('均衡滤波后的I路时域波形');
    % figure(26);
    % [H9_I,w9_I] = freqz(rcsv_msg_source_MF_I, 1, 'whole');
    % plot(w9_I/pi,abs(H9_I));
    % title('均衡滤波后的I路频域波形');

    % figure(27)
    % % plot(rcsv_msg_source_MF_Q);
    % stem(rcsv_msg_source_MF_Q,'.');
    % title('均衡滤波后的Q路时域波形');
    % figure(28);
    % [H9_Q,w9_Q] = freqz(rcsv_msg_source_MF_Q, 1, 'whole');
    % plot(w9_Q/pi,abs(H9_Q));
    % title('均衡滤波后的I路频域波形');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 抽样 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decision_site = fir_M / 2 + span * sps + 1;
    % display(decision_site);
    rcsv_msg_source_MF_option_I = rcsv_msg_source_MF_I(decision_site:sps:end);
    rcsv_msg_source_MF_option_Q = rcsv_msg_source_MF_Q(decision_site:sps:end);
    % figure(29)
    % % plot(rcos_msg_source_addnoise);
    % stem(rcsv_msg_source_MF_option_I,'.');
    % title('I路经过延时处理并抽取后的接收时域波形');
    % figure(30);
    % [H10_I,w10_I] = freqz(rcsv_msg_source_MF_option_I, 1, 'whole');
    % plot(w10_I/pi,abs(H10_I));
    % title('I路经过延时处理并抽取后的接收频域波形');
    % figure(31)
    % % plot(rcos_msg_source_addnoise);
    % stem(rcsv_msg_source_MF_option_Q,'.');
    % title('Q路经过延时处理并抽取后的接收时域波形');
    % figure(32);
    % [H10_Q,w10_Q] = freqz(rcsv_msg_source_MF_option_Q, 1, 'whole');
    % plot(w10_Q/pi,abs(H10_Q));
    % title('Q路经过延时处理并抽取后的接收频域波形');
    msg_source_MF_option_sign_I = sign(rcsv_msg_source_MF_option_I);
    msg_source_MF_option_sign_Q = sign(rcsv_msg_source_MF_option_Q);
    % figure(33);
    % stem(msg_source_MF_option_sign_I, '-*');
    % title('I路判决结果');
    % figure(34);
    % stem(msg_source_MF_option_sign_Q, '-*');
    % title('Q路判决结果');

    %%%%%%%%%%%%%%%%%%% 合并 I路 和 Q路
    length_I = length(msg_source_MF_option_sign_I);
    length_Q = length(msg_source_MF_option_sign_Q);
    % 初始化结果数组
    bipolar_combinedArray = zeros(1, length_I + length_Q);
    % 交替合并
    for i = 1:max(length_I, length_Q)
        if i <= length_I
            bipolar_combinedArray(2*i - 1) = msg_source_MF_option_sign_I(i);
        end
        if i <= length_Q
            bipolar_combinedArray(2*i) = msg_source_MF_option_sign_Q(i);
        end
    end
    combinedArray = (bipolar_combinedArray + 1) / 2;
    % 结果数组
    % disp(combinedArray);
    % stem(combinedArray, '-*');

    %%%%%%%%%%%%%%%%%%%%%%%%% 信宿 %%%%%%%%%%%%%%%%%%%%%%%%%%%
    [err_number(j), bit_err_ratio(j)]=biterr(msg_source(1:length(combinedArray)), combinedArray);
end

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 仿真结果 %%%%%%%%%%%%%%%%%%%%%%%%%
ber_bpsk = berawgn(ebn0, 'psk', 2, 'nondiff');
ber_qpsk = berawgn(ebn0, 'psk', 4, 'nondiff');
semilogy(ebn0, bit_err_ratio, '-*', ebn0, ber_bpsk, '-+', ebn0, ber_qpsk, '-o');
xlabel('比特信噪比');
ylabel('误码率');
title('不同信噪比下误码率仿真曲线');
legend('实验QPSK曲线 ', '理论BPSK曲线', '理论QPSK曲线')
grid on;

