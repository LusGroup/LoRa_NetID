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
sig = phy.modulate(symbols);
% 在signal前后加一些底噪
% pre = zeros(round(0.2*phy.sample_num_per_symbol),1);
% pre = awgn(pre,10,0);
% sig = [pre;sig;pre];

% save sig to file 
filename =  "lora_packet.cfile";
LoRaPHY.write(sig,filename);

% read sig from file
sig = LoRaPHY.read_file(filename);
LoRaPHY.plot_timefrequency(sig,fs,sf,bw);
[symbols_d, cfo, netid] = phy.demodulate(sig);
[data, checksum] = phy.decode(symbols_d);
fprintf("[decode] data:\n");
LoRaPHY.print_payload(data);