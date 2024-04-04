%% 本仿真基于IEEE802.11a
clc;
clear;
close all;
%% 原始序列输入 
num_in = round(rand(1,1e4));  
%% 参数设置
%调制参数设置
M = 4;      %调制阶数。1-BPSK，2-QPSK,3-8QAM,4-16QAM，5-32QAM,6-64QAM
k = 40 ;    %OFDM符号数
code_rate = 2; %卷积编码速率。6--无卷积编码，2--1/3速率，3--1/2速率，4--2/3速率
leng_num_in = k .* 48 .* M ./6 .* code_rate;
num_in = num_in(1:leng_num_in);     %截取输入序列长度
w = length(num_in)./code_rate.*6 /192;      %调制前数据组数目
%信道参数设置 
SNR = 10;   %信噪比
fc = 70e6;  %载波频率
fs = 200e6;     %采样频率
phase_py = 0;    %载波相偏
freq_py = 0;     %载波频偏
%莱斯信道设置
ricianChan = comm.RicianChannel(...
        "SampleRate",             20e6,...
        'PathDelays',             [1.7*50e-9, 2.8*50e-9], ...
        'AveragePathGains',       [-5,-8], ...
        "KFactor",                10,...
        'NormalizePathGains',     true, ... 
        "DirectPathDopplerShift", 20,...
        "MaximumDopplerShift",    10,...
        "RandomStream",           "mt19937ar with seed",...
        "Seed",                   38);
%% 扰码
scram_seed0 = [1,0,1,1,1,0,1];       %扰码寄存器初值
scramnler = scram_seed0;             %扰码寄存器                                                      %数据个数        
scram_in = num_in;        %产生输入随机序列
scram_out0 = zeros(1,length(num_in));             %初始化输出序列
for m = 1:length(scram_in)
    scram_out0(m) = mod(scramnler(1) + scramnler(4) + scram_in(m), 2);       %扰码：7+4+输入数据
    scramnler(:,1:end) = [scramnler(:,2:end), mod(scramnler(1) + scramnler(4), 2)];     %扰码寄存器移位，最低位为7+4
end
%% 卷积编码
conv_in = scram_out0;    
L = 7;          %卷积编码约束长度
switch(code_rate)   %卷积编码速率控制
    case 6  
        conv_out = conv_in;
    case 2  
        puncpat = [1;1;1;1];    %打孔方式
        trellis = poly2trellis(L,[133,171,165]);
        conv_out = convenc(conv_in,trellis,puncpat);
    case 3
        puncpat = [1;1;1;1];
        trellis = poly2trellis(L,[133,171]);
        conv_out = convenc(conv_in,trellis,puncpat);
    case 4
        puncpat = [1;1;1;0];
        trellis = poly2trellis(L,[133, 171]);
        conv_out = convenc(conv_in,trellis,puncpat);
    otherwise
        disp('code_rate_error');
end
%% 一级交织
%一级交织器产生
list = 12;
row = length(conv_out)/list;             %将数据转为12列的矩阵                      
ram = zeros(row,list);              %将输入数据存储在12列，row行的矩阵中
for n = 1:w        %将ram矩阵拆为n个16*12的矩阵，每次写入16行
    for m = 1:192                   %以行写入
        row_index = ceil(m/12);         %写入到矩阵中行的位置
        list_index = mod((m-1),12) + 1;     %写入到矩阵中列的位置
        ram((n-1)*16+row_index,list_index) = conv_out((n-1)*192+m);  %将数据写入ram矩阵
    end
end
int_lea_1_out = zeros(1,length(conv_out));
for n = 1:w
    for list_index = 1:12
        for row_index = 1:16        %按列读出ram中的数据，每次读16行
            int_lea_1_out((n-1)*192+(list_index-1)*16+row_index) = ...
                ram(((n-1)*16+row_index),list_index);
        end
    end
end
%% 二级交织
int_lea_2_out = int_lea_1_out;
for index = 1:192*w         %数据前12个不变，接下来12个两两交换位置
        if(mod((ceil(index/12)-1),2)==1)    %判断数据是不是在后12个位置
            if(mod(index,2) == 0)           %在偶数位置时，将前一个位置的数据给他
                int_lea_2_out(index) = int_lea_1_out(index-1);
            else                            %在奇数位置时，将后一个位置的数据给他
                int_lea_2_out(index) = int_lea_1_out(index+1);
            end
        else
            int_lea_2_out(index) = int_lea_1_out(index);    
        end
end
%% 调制映射
mod_out = qammod(int_lea_2_out', 2^M, 'InputType', 'bit', 'UnitAveragePower', true, 'PlotConstellation', true)';
%% 插入导频到6，20，33，47
%插入导频极性控制，扰码
scram_seed2 = [1,1,1,1,1,1,1];
scram_reg = scram_seed2;
scram_out = zeros(1,k);
for m = 1:k
    scram_out(m) = mod(scram_reg(1) + scram_reg(4), 2);
    scram_reg(:,1:end) = [scram_reg(:,2:end), mod(scram_reg(1) + scram_reg(4), 2)];     %扰码寄存器移位，最低位为7+4
end
rx_interFrq = mod_out;
interFrq_out = zeros(1,52*k);
for m = 1:k
    reg48 = rx_interFrq((m-1)*48+1:m*48);
    reg_pn = scram_out(m);
        if(reg_pn)
            reg_interFrq = [1,1,1,-1];
        else
            reg_interFrq = [-1,-1,-1,1];
        end
        reg52(1:5) = reg48(1:5);
        reg52(7:19) = reg48(6:18);
        reg52(21:32) = reg48(19:30);
        reg52(34:46) = reg48(31:43);
        reg52(48:52) = reg48(44:48);
        reg52(6) = reg_interFrq(1);
        reg52(20) = reg_interFrq(2);
        reg52(33) = reg_interFrq(3);
        reg52(47) = reg_interFrq(4);
        interFrq_out((m-1)*52+1:m*52) = reg52;
end
%% IFFT
ifft_in = interFrq_out;
reg64 = zeros(1,64);
reg_ifft = zeros(1,64);
ifft_out = zeros(1,k*64);
for m = 1:k
    reg52 = ifft_in((m-1)*52+1:m*52);
    reg64(1) = 0;
    reg64(39:64) = reg52(1:26);
    reg64(2:27) = reg52(27:52);
    reg64(28:38) = zeros(1,11);
    reg_ifft = 64^0.5*ifft(reg64);
    ifft_out((m-1)*64+1:m*64) = reg_ifft;
end
%% 加循环前缀CP，加窗
add_cp_in = ifft_out;
add_cp_out = zeros(1,80*k);
reg80 = zeros(1,80);
tmp_end = zeros(1,k+1);
for m = 1:k
    reg64 = add_cp_in((m-1)*64+1:m*64);
    reg80(1:16) = reg64((end-15):end);
    reg80(17:end) = reg64(1:end);
    tmp_end(m+1) = 0.5*reg64(1);
    add_cp_out((m-1)*80+1:m*80) = [0.5*tmp_end(m)+0.5*reg80(1),reg80(2:end)] ;
end
%% 加帧头
%产生长短训练序列、帧头
short_training =[0,0,0,0,-1-1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,1+ ...
    1i,0,0,0,1+1i,0,0,0,1+1i,0,0,0,0,0,0,0,0,0,0,0,0,0,...
    0,0,1+1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,-1-1i, ...
    0,0,0,-1-1i,0,0,0,1+1i,0,0,0]; % short training sequence
sts_frq = (13/6)^0.5 .* short_training;
sts_time = 64^0.5*ifft(sts_frq,64);    %对短训练序列进行ifft
%取16点，将该序列重复10次
sts16 = sts_time(1:16);
sts_rom = sts16;
for n = 1:9
    sts_rom = [sts_rom,sts16];
end
sts_rom = [0.5*sts_rom(1),sts_rom(2:end),0.5*sts16(1)];     %加窗
lts = [ 0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1, ...
    -1,1,-1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,...
    -1,-1,1,1,-1,1,-1,1,1,1,1 ];
lts_time = 64^0.5*ifft(lts,64);        %长训练序列进行ifft
lts_rom = [0.5*lts_time(33),lts_time(34:end),lts_time,lts_time,0.5*lts_time(1)];    %长训练序列加窗
preamble = [sts_rom(1:(end-1)),sts_rom(end)+lts_rom(1),lts_rom(2:end)];     %帧头加窗      
%合成帧头与数据
preamble_out = [preamble(1:(end-1)),preamble(end)+add_cp_out(1),add_cp_out(2:end)];
%% 加帧头噪声
%计算帧头信号均值与方差
mean_preamble = mean(preamble);
std_preamble = std2(preamble);
%产生与帧头能量一致的噪声
noise_pre = mean_preamble + std_preamble*0.707*randn(1,200) + std_preamble*0.707j*randn(1,200);
noise_pre_out = [noise_pre,preamble_out];   %添加帧头噪声
%前后各添加200点0值
zeros200 = zeros(1,200);
pre_zero_out = [zeros200,noise_pre_out,zeros200];
%% 信道。加频偏、相偏
%基带成型滤波，内插10倍
rcosine_filter = rcosdesign(0.3,15,10,'sqrt');
FIRInter_rcosine = dsp.FIRInterpolator('InterpolationFactor',10,'Numerator',rcosine_filter); 
rcosine_filter_out = step(FIRInter_rcosine,pre_zero_out')';
%载波调制
t = 0 : 1/fs : (length(rcosine_filter_out)-1)/fs;
I_carry = cos(2*pi*fc*t);
Q_carry = -sin(2*pi*fc*t);
tx_carry_out = real(rcosine_filter_out).*I_carry + imag(rcosine_filter_out).*Q_carry;   %载波调制
%经过莱斯信道
rician_chan_out = ricianChan(tx_carry_out')';
%经过AWGN信道
awgn_out = awgn(rician_chan_out, SNR, 'measured');
%解载波调制，减采样10倍
I_carry = cos(2*pi*(fc+freq_py)*t+phase_py);
Q_carry = -sin(2*pi*(fc+freq_py)*t+phase_py);

rx_de_carry_out = 2*awgn_out.*I_carry + 2j*awgn_out.*Q_carry;

down_sample_filter = dsp.FIRDecimator('DecimationFactor',10,'Numerator',rcosine_filter);
down_carry_out = step(down_sample_filter ,rx_de_carry_out')';
%% 延时相关运算
rx_relation = down_carry_out;
reg16 = zeros(1,16);
reg16_delay = zeros(1,16);
n_rela = 800;
judge_m = zeros(1,n_rela-16);
for m = 17:n_rela
    reg16_delay = rx_relation(m-16:m-1);
    reg16(1) = (conj(reg16_delay(1)).*rx_relation(m));
    rela = abs(sum(reg16));
    power = sum(abs(reg16_delay).^2);
    judge_m(m-16) = rela/power;
    reg16 = [0,reg16(1:(end-1))];
end
figure,plot(judge_m);
title("帧同步判决变量");
%% 长度保持
rx_len_hold = judge_m;
reg32 = zeros(1,32);
detect_en = zeros(1,length(rx_len_hold)-48);
detect = 0;
for m = 1:length(rx_len_hold)-48
    reg32 = rx_len_hold(m:m+31);
    reg48 = rx_len_hold(m:m+47);
    if(reg32 > 0.5)
        detect = 1;
    end
    if(reg48 < 0.5)
        detect = 0;
    end
    if(detect)
        detect_en(m) = 1;
    else
        detect_en(m) = 0;
    end
end
detect_index = find(detect_en); 
if ~isempty(detect_index)
    detect_index_first = detect_index(1);
else
    disp('帧同步失败');
end
%% 载波同步
%延迟相关，相关累加
reg_sum = zeros(1,4);
for n = 1:4
    reg16_1 = rx_relation(detect_index_first+16*(n-1):detect_index_first+16*n-1);
    reg16_2 = rx_relation(detect_index_first+16*n:detect_index_first+16*(n+1)-1);
    reg16 = conj(reg16_1) .* reg16_2;
    reg_sum(n) = (sum(reg16));
end
%% 频偏计算
angle_sum = zeros(1,length(reg_sum));
for n = 1:length(reg_sum)
    angle_sum(n) = angle(reg_sum(n));
end
ave_angle = sum(angle_sum)./4;
ave_angle = ave_angle./16;
%偏移累加
reg_phase = zeros(1,length(rx_relation));
for n = 2:length(reg_phase)
    reg_phase(n) = reg_phase(n-1) - ave_angle;
    if(reg_phase(n) < -1*pi)
        reg_phase(n) = reg_phase(n) + 2*pi;
    end
    if(reg_phase(n) > 1*pi)
        reg_phase(n) = reg_phase(n) - 2*pi;
    end
end
%补偿因子计算
compensated_data = 1j*sin(reg_phase) + cos(reg_phase);
reg_comp_data = rx_relation((detect_index_first):end); 
reg_comp_out = reg_comp_data .* compensated_data(1:length(reg_comp_data));  %频偏补偿
%% 符号同步
match_energy = zeros(1,160);
for m = 1:160-15
    match_energy(m) = abs(sum(reg_comp_out(m:m+15) .* conj(sts16)));
end
figure,stem(match_energy);
title("符号同步能量");
peak_cnt = 0;
peak_order = sort(match_energy);    %峰值升序排列
judge_energy = peak_order(end-8);   %取第9大峰值为门限
for m = 1:length(match_energy)
    if(match_energy(m) >= judge_energy)
        peak_cnt = peak_cnt + 1;
        if(peak_cnt == 9)   %9个能量尖峰
            lts_start_index = m + 12;   %定位位置调整,变化范围为1-16.CP长度为16。16为无相位偏转精准定位。
            break;
        end
    end
end
% 提取长训练序列
rx_lts = reg_comp_out(lts_start_index : lts_start_index + 160 - 1);
%% 解循环前缀
rx_de_cp = reg_comp_out(lts_start_index + 160 : 80*k + lts_start_index + 160);
for n = 1:k
    reg80 = rx_de_cp((n-1)*80+1:n*80);
    reg64 = reg80(17:80); 
    de_add_cp_out((n-1)*64+1:n*64) = reg64;
end
%% fft
rx_fft = de_add_cp_out;
fft_out = zeros(1,k*64);
reg_fft = zeros(1,64);
for n = 1:k
    reg64 = rx_fft((n-1)*64+1:n*64); 
    reg_fft = 64^-0.5*fft(reg64);
    fft_out((n-1)*64+1:n*64) =  reg_fft;
end
%% 信道估计
rx_lts1 = 64^-0.5*(fft(rx_lts(33 : 33+64-1)));
rx_lts2 = 64^-0.5*(fft(rx_lts(33+64 : end)));
R_rlts = (rx_lts1 + rx_lts2) ./ 2;
H = R_rlts .* conj(lts);
E_rtl = R_rlts .* conj(R_rlts);
H_out = zeros(1,length(fft_out));
for n = 1:k
    H_out((n-1)*64+1 : n*64) =  fft_out((n-1)*64+1 : n*64) .* conj(H) ./E_rtl;
end
%% 去空子载波
rx_de_sub = H_out;
de_sub_out = zeros(1,k*52);
for n = 1:k
    reg64 = rx_de_sub((n-1)*64+1:n*64); 
    reg52(1:26) = reg64(39:64);
    reg52(27:52) = reg64(2:27);
    de_sub_out((n-1)*52+1:n*52) =  reg52;
end
%% 解插入导频
rx_interFrq = de_sub_out;
de_interFrq_out = zeros(1,48*k);
rx_reg_interFrq = zeros(1,4*k);
for n = 1:k
    reg52 = rx_interFrq((n-1)*52+1:n*52);
    reg48(1:5) = reg52(1:5);
    reg48(6:18) = reg52(7:19);
    reg48(19:30) = reg52(21:32);
    reg48(31:43) = reg52(34:46);
    reg48(44:48) = reg52(48:52);
    de_interFrq_out((n-1)*48+1:n*48) = reg48;
    rx_reg_interFrq(1+(n-1)*4) = reg52(6);
    rx_reg_interFrq(2+(n-1)*4) = reg52(20);
    rx_reg_interFrq(3+(n-1)*4) = reg52(33);
    rx_reg_interFrq(4+(n-1)*4) = reg52(47);
end
figure,scatter (real(de_interFrq_out),imag(de_interFrq_out),'filled') %星座图绘制 
title('剩余相位追踪前');
%% 剩余相位追踪
sfo_comp_out = zeros(1,length(de_interFrq_out));
for n = 1:k
    reg_pn = scram_out(n);
    if(reg_pn)
        reg_interFrq = [1,1,1,-1];
    else
        reg_interFrq = [-1,-1,-1,1];
    end
    phase_factor = conj(sum(reg_interFrq .* rx_reg_interFrq((n-1)*4+1 : n*4))) ./ 4;    %相位因子计算
    sfo_comp_out((n-1)*48+1:n*48) = de_interFrq_out((n-1)*48+1:n*48) .* phase_factor;   %相位补偿
end
%% qam解调
rx_demod = sfo_comp_out;
scatterplot(sfo_comp_out);
title('剩余相位追踪后');
demod_out = qamdemod(rx_demod',2^M,'UnitAveragePower', true,'OutputType','bit')';
%% 解二级交织
de_int_lea_2_out = demod_out;       %与二级交织过程相同
for index = 1:length(de_int_lea_2_out)
        if(mod((ceil(index/12)-1),2)==1)
            if(mod(index,2) == 0)
                de_int_lea_2_out(index) = demod_out(index-1);
            else
                de_int_lea_2_out(index) = demod_out(index+1);
            end
        else
            de_int_lea_2_out(index) = demod_out(index);
        end   
end

%% 解一级交织
rx_int_lea_1 = zeros(list,row);      %将收到的数据写入12行，row列的矩阵中，与一级交织时相反
for n = 1:w         %将数据以192为单位，分割为n份
    for m = 1:192
        row_index = ceil(m/16);      %写入的行的位置
        list_index = mod((m-1),16) + 1;      %写入的列的位置
        rx_int_lea_1(row_index,list_index+(n-1)*16) = de_int_lea_2_out((n-1)*192+m);     %按列写入，总共12行，16列
    end
end
de_int_lea_1_out = zeros(1,length(de_int_lea_2_out));
for n = 1:w         %将矩阵中的数据按行读出
    for list_index = 1:16
        for row_index = 1:12        %按行读出矩阵数据，一共12行
            de_int_lea_1_out((n-1)*192+(list_index-1)*12 + row_index) = ...
                rx_int_lea_1(row_index,(n-1)*16+list_index);         %每读出16列后，数据位数加上192
        end
    end
end
%% 解卷积编码
L = 7;
tblen = 6*L;
switch(code_rate)
    case 6
        vtb_out = de_int_lea_1_out;
    case 2
        puncpat = [1;1;1;1];
        trellis = poly2trellis(L,[133,171,165]);
        vtb_out = vitdec(de_int_lea_1_out,trellis,tblen,'trunc','hard',puncpat);
    case 3
        puncpat = [1;1;1;1];
        trellis = poly2trellis(L,[133,171]);
        vtb_out = vitdec(de_int_lea_1_out,trellis,tblen,'trunc','hard',puncpat);
    case 4
        puncpat = [1;1;1;0];
        trellis = poly2trellis(L,[133, 171]);
        vtb_out = vitdec(de_int_lea_1_out,trellis,tblen,'trunc','hard',puncpat);
    otherwise
        disp('code_rate_error');
end
%% 解扰码
scramnler_rx = scram_seed0;      %解扰码与扰码方法一致
de_scram_out0 = zeros(1,length(vtb_out));
for n = 1:length(vtb_out)
    de_scram_out0(n) = mod(scramnler_rx(1,1) + scramnler_rx(1,4) + vtb_out(n), 2); 
    scramnler_rx(:,1:end) = [scramnler_rx(:,2:end), mod(scramnler_rx(1) + scramnler_rx(4), 2)];
end
%% 差错
N_error = 0;
for n = 1:length(de_scram_out0)
    if(de_scram_out0(n) ~= num_in(n))
        N_error = N_error + 1;
    end    
end 
bit_error_rate = N_error/length(de_scram_out0)%误码率结果