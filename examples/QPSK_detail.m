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
% msg_source = [ones(1,20) zeros(1,20) randi([0,1], 1, 99960)];
% msg_source = [ones(1,10) zeros(1,10) randi([0,1], 1, 10)];

msg_source = [1 0 1 0 1 0 1 0 1 0 1 0 1 1 1 1 1 0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%   发射机  %%%%%%%%%%%%%%%%%%%%
bipolar_msg_source = 2*msg_source-1;

figure(1);
stem(bipolar_msg_source, '.');axis([0,length(msg_source),-2,2]);
title('信源时域波形');

bipolar_msg_source_I = bipolar_msg_source(1:2:end);
bipolar_msg_source_Q = bipolar_msg_source(2:2:end);

% display(bipolar_msg_source_I);
% display(bipolar_msg_source_Q);





