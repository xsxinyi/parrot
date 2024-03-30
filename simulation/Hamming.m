clc
clear
close all
% Title:  (7,4)汉明码在BPSK调制系统中的性能 %
% Data:    2020.01.20 %
% Author:  K.X.Song   % 
L_data = 1000000;                       % 原始数据长度
L_code = L_data/4*7;                    % 编码后数据长度
data = round(rand(1,L_data));           % 原始数据
EbN0_dB = 0:10;                         % Eb/N0 dB形式
EbN0 = 10.^(EbN0_dB/10);                % Eb/N0
Eb = 7/4;                               % 每比特能量
N0 = Eb ./ EbN0;                        % 噪声功率
error = zeros(1,length(EbN0_dB));       % 预置错误个数
ber = zeros(1,length(EbN0_dB));         % 预置仿真误比特率
tber_BPSK = zeros(1,length(EbN0_dB));   % 预置BPSK理论误比特率
tber_BPSK_Hamming = zeros(1,length(EbN0_dB));   % 预置BPSK+(7,4)汉明码理论误比特率
tber = zeros(1,length(EbN0_dB));       % 预置理论误比特率
% (7,4)汉明码
code = zeros(1,L_code);                  % 预置编码后的序列
for k = 0:L_data/4-1
    code(1,k*7+1) = data(1,k*4+1);      % 编码方法详见技术文档
    code(1,k*7+2) = data(1,k*4+2);
    code(1,k*7+3) = data(1,k*4+3);
    code(1,k*7+4) = data(1,k*4+4);
    code(1,k*7+5) = xor(data(1,k*4+2),xor(data(1,k*4+3),data(1,k*4+4)));    
    code(1,k*7+6) = xor(data(1,k*4+1),xor(data(1,k*4+3),data(1,k*4+4)));
    code(1,k*7+7) = xor(data(1,k*4+1),xor(data(1,k*4+2),data(1,k*4+4)));
end
send = (code - 1/2) * 2;                % BPSK调制
for q = 1:length(EbN0_dB)
    noise = sqrt(N0(q)/2) * randn(1,L_code);    % AWGN
    receive = send + noise;                     % 接收数据
    demod = zeros(1,L_code);                    % 预置解调数据
    decode = zeros(1,L_data);                   % 预置解码数据
    for w = 1:L_code
        if (receive(w) >= 0)
            demod(w) = 1;              % 数轴右侧 ->  1
        else
            demod(w) = 0;              % 数轴左侧 ->  0
        end
    end
    for k = 0:L_data/4-1
        S(1) = xor(demod(1,k*7+2),xor(demod(1,k*7+3),xor(demod(1,k*7+4),demod(1,k*7+5))));    % 校正子详见技术文档
        S(2) = xor(demod(1,k*7+1),xor(demod(1,k*7+3),xor(demod(1,k*7+4),demod(1,k*7+6))));
        S(3) = xor(demod(1,k*7+1),xor(demod(1,k*7+2),xor(demod(1,k*7+4),demod(1,k*7+7))));
        if (S == [0,1,1])
            demod(1,k*7+1) = ~demod(1,k*7+1);
        elseif (S == [1,0,1])
            demod(1,k*7+2) = ~demod(1,k*7+2);
        elseif (S == [1,1,0])
            demod(1,k*7+3) = ~demod(1,k*7+3);
        elseif (S == [1,1,1])
            demod(1,k*7+4) = ~demod(1,k*7+4);
        elseif (S == [1,0,0])
            demod(1,k*7+5) = ~demod(1,k*7+5);
        elseif (S == [0,1,0])
            demod(1,k*7+6) = ~demod(1,k*7+6);
        elseif (S == [0,0,1])
            demod(1,k*7+7) = ~demod(1,k*7+7);
        end
        decode(1,k*4+1) = demod(1,k*7+1);
        decode(1,k*4+2) = demod(1,k*7+2);
        decode(1,k*4+3) = demod(1,k*7+3);
        decode(1,k*4+4) = demod(1,k*7+4);
    end
    for w = 1:L_data
        if (decode(w) ~= data(w))
            error(q) = error(q) + 1;    % 错误比特个数
        end
    end
    ber(q) = error(q) / L_data;         % 仿真误比特率
    % 理论误比特率的推导见: 《通信原理》第47页
    % 对于实数, delta^2 = n0 / 2, 而 A = sqrt(Eb) 的. 所以 A/delta =
    % sqrt(Eb)/sqrt(n0/2) = sqrt(2*Ebn0)
    tber(q) = qfunc(sqrt(2*EbN0(q)));   % 理论误比特率
end
figure
semilogy(EbN0_dB,ber,'M-X',EbN0_dB,tber,'B-O');     % 画图
grid on;                                        % 坐标轴开启
axis([0 10 10^-5 10^-1])                        % 限制作图范围
xlabel('Eb/N0 (dB)');                           % 横坐标
ylabel('BER');                                  % 纵坐标
legend('BPSK+(7,4)汉明码仿真误比特率','BPSK理论误比特率');   % 图例
