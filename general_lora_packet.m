close all;
clear all;
freq = 868e6;
bw = 250e3;
fs = 500e3;
sf = 7;       % sampling rate 1 MHz

phy = LoRaPHY(freq,sf,bw,fs);
phy.has_header = 1;
phy.CR = 3;
phy.CRC = 1;
phy.preamble_len = 6;  
% preamble: 8 basic upchirps
% data = phy.Ref_upchirp;

% LoRa_Analyse.plot_timefrequency(phy.Ref_upchirp,phy.Fs,phy.SF,phy.BW);
% data is 'hello world,this is a message by simulating in maltab'
data = [104 101 108 108 111 32 119 111 114 108 100 44 116 104 105 115 32 105 115 32 97 32 109 101 115 115 97 103 101 32 98 121 32 115 105 109 117 108 97 116 105 110 103 32 105 110 32 109 97 108 116 97 98];
symbols = phy.encode(data.');

uc = LoRaPHY.chirp(true,sf,bw,fs,0,0);
dc = LoRaPHY.chirp(false,sf,bw,fs,0,0);

preamble = [uc;dc;uc;dc;uc;dc;uc;dc];
% netid = [LoRaPHY.chirp(true,sf,bw,fs,8,0);LoRaPHY.chirp(true,sf,bw,fs,16,0);];
chirp_len = length(uc);
sfd = [dc; dc; dc(1:round(chirp_len/4))];

data = zeros(length(symbols)*chirp_len, 1); 
for i = 1:length(symbols)
    data((i-1)*chirp_len+1:i*chirp_len) =  LoRaPHY.chirp(true, sf, bw, fs, symbols(i), 0); %data数组只能1开始
end
sig = [preamble;sfd; data];
% LoRa_Analyse.plot_timefrequency(sig,phy.Fs,phy.SF,phy.BW);
% 在signal前后加一些底噪
% pre = zeros(5.25*phy.sample_num_per_symbol,1);
% pre = awgn(pre,10,0);
% sig = [pre;sig;pre];
filename =  "../sig/ud_pre_packet.cfile";
LoRa_Analyse.write(sig,filename);