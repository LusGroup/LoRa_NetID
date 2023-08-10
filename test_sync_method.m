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
filename =  "../sig/ud_pre_packet.cfile";
sig = LoRaPHY.read_file(filename);
LoRa_Analyse.plot_timefrequency(sig,fs,sf,bw);
[symbols_d, cfo, netid] = phy.demodulate(sig);
[data, checksum] = phy.decode(symbols_d);
fprintf("[decode] data:\n");
LoRaPHY.print_payload(data);