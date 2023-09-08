classdef LoRaPHY  < handle & matlab.mixin.Copyable

    properties
        C_freq                     % 中心频率
        SF                         % 扩频因子（7,8,9,10,11,12）        
        BW                         % 带宽
        Fs                         % 采样频率 
        CR                         % 码率：（1:4/5 2:4/6 3:4/7 4:4/8）
        Ts                         % LoRa符号持续时间 2^(sf)/bw
        Payload_len                % 有效载荷的长度
        preamble_len               % 前导部分的长度
        CFO                        % 载波频率偏移
        has_header                 % 隐式header：0 显示header：1
        CRC                        % CRC = 1 如果开启CRC检验
        LDR                        % LDR = 1 如果开启 low data rate 模式
        
        Whitening_seq              % 白化序列
        Crc_generator              % CRC 产生器
        Header_checksum_matrix     % header CRC校验矩阵
        is_debug                   % is_debug = true 如果需要显示debug信息
        hamming_decoding_en        % hamming_decoding_en = 1 如果开启hanmming解码
        
        fft_len                    % fft size
        bin_num                    % fft后bin的个数（包括Zero-padding）
        preamble_bin               % 当前解调窗口下前导部分的bin，用来消除cfo的影响
        sample_num_per_symbol      % 一个符号的采样点个数
        zero_padding_ratio         % fft零比特填充的比例
        
        fast_mode                  % set `true` for fast execution (ignore low-pass filter)
        
        Ref_upchirp                % 参考base upchirp
        Ref_downchirp              % 参考base downchirp
        Ref_pk_bin_list            % 参考fft bin list
        
        sig                        % 输入的基带信号
        
        m_n_up_req
        os_factor                  % fs/bw的比值

        pattern                    % 前导码的形式
        pattern_bin_list           % 前导码对应的fft bin
        
    end

    methods
     
        function self = LoRaPHY(C_freq,SF,BW,Fs,pattern)

            if nargin < 5
                self.pattern = [];
            else
                self.pattern = pattern;
            end
        
            self.C_freq = C_freq;
            self.SF = SF;
            self.Fs = Fs;
            self.BW = BW;
            self.Ts = 2^(SF)/BW;
            self.zero_padding_ratio = 10;
            
            self.has_header = 1;
            self.CR = 1;
            self.CRC = 1;
            self.is_debug = false;
            self.hamming_decoding_en = true;
            self.fast_mode = false;
            self.CFO = 0;
            self.pattern_bin_list = [];
            
            self.Whitening_seq = uint8([
                0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE1, 0xC2, 0x85, 0x0B, 0x17, 0x2F, 0x5E, 0xBC, 0x78, 0xF1, 0xE3, ...
                0xC6, 0x8D, 0x1A, 0x34, 0x68, 0xD0, 0xA0, 0x40, 0x80, 0x01, 0x02, 0x04, 0x08, 0x11, 0x23, 0x47, ...
                0x8E, 0x1C, 0x38, 0x71, 0xE2, 0xC4, 0x89, 0x12, 0x25, 0x4B, 0x97, 0x2E, 0x5C, 0xB8, 0x70, 0xE0, ...
                0xC0, 0x81, 0x03, 0x06, 0x0C, 0x19, 0x32, 0x64, 0xC9, 0x92, 0x24, 0x49, 0x93, 0x26, 0x4D, 0x9B, ...
                0x37, 0x6E, 0xDC, 0xB9, 0x72, 0xE4, 0xC8, 0x90, 0x20, 0x41, 0x82, 0x05, 0x0A, 0x15, 0x2B, 0x56, ...
                0xAD, 0x5B, 0xB6, 0x6D, 0xDA, 0xB5, 0x6B, 0xD6, 0xAC, 0x59, 0xB2, 0x65, 0xCB, 0x96, 0x2C, 0x58, ...
                0xB0, 0x61, 0xC3, 0x87, 0x0F, 0x1F, 0x3E, 0x7D, 0xFB, 0xF6, 0xED, 0xDB, 0xB7, 0x6F, 0xDE, 0xBD, ...
                0x7A, 0xF5, 0xEB, 0xD7, 0xAE, 0x5D, 0xBA, 0x74, 0xE8, 0xD1, 0xA2, 0x44, 0x88, 0x10, 0x21, 0x43, ...
                0x86, 0x0D, 0x1B, 0x36, 0x6C, 0xD8, 0xB1, 0x63, 0xC7, 0x8F, 0x1E, 0x3C, 0x79, 0xF3, 0xE7, 0xCE, ...
                0x9C, 0x39, 0x73, 0xE6, 0xCC, 0x98, 0x31, 0x62, 0xC5, 0x8B, 0x16, 0x2D, 0x5A, 0xB4, 0x69, 0xD2, ...
                0xA4, 0x48, 0x91, 0x22, 0x45, 0x8A, 0x14, 0x29, 0x52, 0xA5, 0x4A, 0x95, 0x2A, 0x54, 0xA9, 0x53, ...
                0xA7, 0x4E, 0x9D, 0x3B, 0x77, 0xEE, 0xDD, 0xBB, 0x76, 0xEC, 0xD9, 0xB3, 0x67, 0xCF, 0x9E, 0x3D, ...
                0x7B, 0xF7, 0xEF, 0xDF, 0xBF, 0x7E, 0xFD, 0xFA, 0xF4, 0xE9, 0xD3, 0xA6, 0x4C, 0x99, 0x33, 0x66, ...
                0xCD, 0x9A, 0x35, 0x6A, 0xD4, 0xA8, 0x51, 0xA3, 0x46, 0x8C, 0x18, 0x30, 0x60, 0xC1, 0x83, 0x07, ...
                0x0E, 0x1D, 0x3A, 0x75, 0xEA, 0xD5, 0xAA, 0x55, 0xAB, 0x57, 0xAF, 0x5F, 0xBE, 0x7C, 0xF9, 0xF2, ...
                0xE5, 0xCA, 0x94, 0x28, 0x50, 0xA1, 0x42, 0x84, 0x09, 0x13, 0x27, 0x4F, 0x9F, 0x3F, 0x7F]');
            
            self.Header_checksum_matrix = gf([
            1 1 1 1 0 0 0 0 0 0 0 0
            1 0 0 0 1 1 1 0 0 0 0 1
            0 1 0 0 1 0 0 1 1 0 1 0
            0 0 1 0 0 1 0 1 0 1 1 1
            0 0 0 1 0 0 1 0 1 1 1 1
            ]);
            
            self.Crc_generator = comm.CRCGenerator('Polynomial','X^16 + X^12 + X^5 + 1');
            self.preamble_len = 6;
            self.init();
        end

        function init(self)
            % 初始化参数
            self.os_factor = self.Fs/self.BW;
            self.bin_num = 2^self.SF * self.zero_padding_ratio;
            self.sample_num_per_symbol = 2*2^self.SF;
            self.fft_len = self.sample_num_per_symbol * self.zero_padding_ratio;
            
            self.Ref_upchirp = LoRaPHY.chirp(true,self.SF,self.BW,2*self.BW,0,self.CFO,0);
            self.Ref_downchirp = LoRaPHY.chirp(false,self.SF,self.BW,2*self.BW,0,self.CFO,0);
            self.Ref_pk_bin_list = [];
            fft_resolution = 1/(2*self.Ts)/self.zero_padding_ratio;
            up_fft_bin = self.BW/2/fft_resolution;
            down_fft_bin = 3*self.BW/2/fft_resolution;
            for i = 1:length(self.pattern)
                if (self.pattern(i) == 1)
                    self.Ref_pk_bin_list = [self.Ref_pk_bin_list up_fft_bin];
                else
                    self.Ref_pk_bin_list = [self.Ref_pk_bin_list down_fft_bin];
                end
            end
            % Low Data Rate Optimization (LDRO) mode in LoRa
            % If the chirp peird is larger than 16ms, the least significant
            % two bits are considered unreliable and are neglected.
            if 2^(self.SF)/self.BW > 16e-3
             self.LDR = 1;
            else
             self.LDR = 0;
            end
        end

        function symbols = encode(self,payload)
            % 该函数用来将输入的 payload 编码成 symbols
            
            if self.CRC
             data = uint8([payload;self.calc_crc(payload)]);
            else
             data = uint8(payload);
            end
            
            plen = length(payload); %只要 payload
            sym_num = self.calc_sym_num(plen);
            
            nibble_num = self.SF - 2 + (sym_num-8)/(self.CR+4)*(self.SF-2*self.LDR);
            
            data_w = uint8([data; 255*zeros(ceil((nibble_num-2*length(data))/2), 1)]); % 多余的nibble用于交织，初始化为0
            
            data_w(1:plen) = self.whitening(data_w(1:plen));
            
            
            data_nibbles = uint8(zeros(nibble_num, 1));
            for i = 1:nibble_num
             idx = ceil(i/2);
             if mod(i, 2) == 1
                 data_nibbles(i) = bitand(data_w(idx), 0xf);
             else
                 data_nibbles(i) = bitshift(data_w(idx), -4);
             end
            end
            
            if self.has_header
             header_nibbles = self.gen_header(plen);
            else
             header_nibbles =[];
            end
            
            codewords = self.hamming_encode([header_nibbles; data_nibbles]);
            
            % 交织
            symbols_i = self.diag_interleave(codewords(1:self.SF-2), 8);
            ppm = self.SF - 2*self.LDR;
            rdd = self.CR + 4;
            for i = self.SF-1:ppm:length(codewords)-ppm+1
             symbols_i = [symbols_i; self.diag_interleave(codewords(i:i+ppm-1), rdd)];
            end
            
            symbols = self.gray_decoding(symbols_i);
        end


        function whitened_payload = whitening(self,data)
            % 该函数用来做白化的操作
            len = length(data);
            whitened_payload = bitxor(data(1:len), self.Whitening_seq(1:len));
            self.print_bin("whiten",whitened_payload);
        end

        function header_nibbles = gen_header(self, plen)
             % gen_header  Generate a valid LoRa header
             %
             % input:
             %     plen: Payload length
             % output:
             %     header_nibbles: Header in nibbles
             header_nibbles = zeros(5,1);
             header_nibbles(1) = bitshift(plen, -4);
             header_nibbles(2) = bitand(plen, 15);
             header_nibbles(3) = bitor(2*self.CR,self.CRC);
             header_checksum = self.Header_checksum_matrix * gf(reshape(de2bi(header_nibbles(1:3), 4, 'left-msb')', [], 1));
             x = header_checksum.x;
             header_nibbles(4) = x(1);
             for i = 1:4
                header_nibbles(5) = bitor(header_nibbles(5), x(i+1)*2^(4-i));
             end
        end 
         
        function codewords = hamming_encode(self,nibbles)
            % hamming_encode  Hamming encoding process in LoRa
            %
            % input:
            %     nibbles: Data in nibbles
            % output:
            %     codewords: Data after hamming encoding
            
            nibble_num = length(nibbles);
            codewords = uint8(zeros(nibble_num,1));
            for i = 1:nibble_num
                nibble = nibbles(i);
            
                p1 = LoRaPHY.bit_reduce(@bitxor, nibble, [1 3 4]);
                p2 = LoRaPHY.bit_reduce(@bitxor, nibble, [1 2 4]);
                p3 = LoRaPHY.bit_reduce(@bitxor, nibble, [1 2 3]);
                p4 = LoRaPHY.bit_reduce(@bitxor, nibble, [1 2 3 4]);
                p5 = LoRaPHY.bit_reduce(@bitxor, nibble, [2 3 4]);
                if i <= self.SF - 2
                    % the first SF-2 nibbles use CR=4/8
                    cr_now = 4;
                else
                    cr_now = self.CR;
                end
                switch cr_now
                    case 1
                        codewords(i) = bitor(uint8(p4)*16, nibble);
                    case 2
                        codewords(i) = LoRaPHY.word_reduce(@bitor, [uint8(p5)*32 uint8(p3)*16 nibble]);
                    case 3
                        codewords(i) = LoRaPHY.word_reduce(@bitor, [uint8(p2)*64 uint8(p5)*32 uint8(p3)*16 nibble]);
                    case 4
                        codewords(i) = LoRaPHY.word_reduce(@bitor, [uint8(p1)*128 uint8(p2)*64 uint8(p5)*32 uint8(p3)*16 nibble]);
                    otherwise
                        % THIS CASE SHOULD NOT HAPPEN
                        error('Invalid Code Rate!');
                end
            end
        end

        function symbols_i = diag_interleave(self, codewords, rdd)
            % diag_interleave  Diagonal interleaving
            %
            % input:
            %     codewords: Data in nibbles
            %     rdd: Bits with redundancy
            %          For example, code rate 4/5 means rdd = 5
            % output:
            %     symbols_i: Symbols after diagonal interleaving
            
            tmp = de2bi(codewords, rdd, 'right-msb');
            symbols_i = uint16(bi2de(cell2mat(arrayfun(@(x) circshift(tmp(:,x), 1-x), 1:rdd, 'un', 0))'));
            self.print_bin("Interleave", symbols_i);
        end

        function symbols = gray_decoding(self, symbols_i)
            % gray_decoding  Gray decoding
            %                `gray_decoding` is used in the ENCODING process
            %
            % input:
            %     symbols_i: Interleaved symbols
            % output:
            %     symbols: Final symbols to be modulated in a packet
            
            symbols = zeros(length(symbols_i), 1);
            for i = 1:length(symbols_i)
                num = uint16(symbols_i(i));
                mask = bitshift(num, -1);
                while mask ~= 0
                    num = bitxor(num, mask);
                    mask = bitshift(mask, -1);
                end
                if i <= 8 || self.LDR
                    symbols(i) = mod(num * 4 + 1, 2^self.SF);
                else
                    symbols(i) = mod(num + 1, 2^self.SF);
                end
            end
        end

        function sym_num = calc_sym_num(self, plen)
            % calc_sym_num  Calculate number of symbols
            %
            % input:
            %     plen: Payload length
            % output:
            %     sym_num: Number of symbols
            
            sym_num = double(8 + max((4+self.CR)*ceil(double((2*plen-self.SF+7+4*self.CRC-5*(1-self.has_header)))/double(self.SF-2*self.LDR)), 0));
        end
         
        function checksum = calc_crc(self, data)
            % calc_crc  Calculate payload CRC
            %
            % input:
            %     data: Data in bytes
            % output:
            %     checksum: CRC result
            
            switch length(data)
                case 0
                    checksum = [0; 0];
                case 1
                    checksum = [data(end); 0];
                case 2
                    checksum = [data(end); data(end-1)];
                otherwise
                    input = data(1:end-2);
                    seq = self.Crc_generator(reshape(logical(de2bi(input, 8, 'left-msb'))', [], 1));
                    checksum_b1 = bitxor(bi2de(seq(end-7:end)', 'left-msb'), data(end));
                    checksum_b2 = bitxor(bi2de(seq(end-15:end-8)', 'left-msb'), data(end-1));
                    checksum = [checksum_b1; checksum_b2];
            end
        end

        %------------------------------------------------------------------------
        %
        % 以下是 modulation 的代码
        %
        % ------------------------------------------------------------------------
            
        function s = modulate(self, symbols)
            % modulate  Modulate a baseband signal
            %
            % input:
            %     symbols: A vector of chirp symbols to be modulated
            %              valid symbol range: 0 to 2^sf-1
            % output:
            %     s: A valid LoRa baseband signal
            
            uc = LoRaPHY.chirp(true,self.SF,self.BW,self.Fs,0,self.CFO);
            dc = LoRaPHY.chirp(false,self.SF,self.BW,self.Fs,0,self.CFO);
            preamble = repmat(uc, self.preamble_len, 1);
            netid = [LoRaPHY.chirp(true,self.SF,self.BW,self.Fs,8,self.CFO); LoRaPHY.chirp(true,self.SF,self.BW,self.Fs,16,self.CFO)];
            
            chirp_len = length(uc);
            sfd = [dc; dc; dc(1:round(chirp_len/4))];
            
            data = zeros(length(symbols)*chirp_len, 1); 
            for i = 1:length(symbols)
                data((i-1)*chirp_len+1:i*chirp_len) =  LoRaPHY.chirp(true, self.SF, self.BW, self.Fs, symbols(i), self.CFO); %data数组只能1开始
            end
            s = [preamble; netid; sfd; data];
        end

        %------------------------------------------------------------------------
        %
        % 以下是 demodulation 的代码
        %
        % ------------------------------------------------------------------------
        function [symbols_m, cfo_m, netid_m] = demodulate(self, sig)
            % demodulate  LoRa packet demodulation
            %
            % input:
            %     sig: Baseband signal in complex
            % output:
            %     symbols_m: A matrix containing the demodulated results.
            %                Each column vector represents the symbols of
            %                a successfully demodulated packet.
            %     cfo_m: A vector containing the carrier frequency offset
            %            results. Each element represents the CFO of the
            %            packet in symbols_m.
            
            self.CFO = 0;
            self.init();
            
            if ~self.fast_mode
             sig = lowpass(sig,self.BW/2,self.Fs);
            end
            
            % resample signal with 2*bandwidth
            self.sig = resample(sig, 2*self.BW, self.Fs);
            
            % 画出重采样后的时频图
%             LoRaPHY.plot_timefrequency(self.sig,self.Fs,self.SF,self.BW);
            symbols_m = [];
            cfo_m = [];
            netid_m = [];
            x = 1;
            while x < length(self.sig)
                x = self.detect(x);  
                % x = self.pre_detect(x);
                % 画出 detect 到的点的时频图
%                 LoRaPHY.plot_timefrequency(self.sig(x:x+10*self.sample_num_per_symbol),self.Fs,self.SF,self.BW);
                if x < 0
                    break;
                end

                % 用于调试此时解调窗口位于解码的哪一部分
                if false
                    signal = self.sig(x:x+self.sample_num_per_symbol*2-1);
                    LoRaPHY.plot_timefrequency(signal,self.Fs,self.SF,self.BW);
                end
                
                % align symbols with SFD
                x = self.sync(x);
                % x = self.pre_sync(x);
                %NetId
                pk_netid1 = self.dechirp(round(x - 4.25 * self.sample_num_per_symbol));
                pk_netid2 = self.dechirp(round(x - 3.25 * self.sample_num_per_symbol));
                netid_m = [netid_m;
                  [mod((pk_netid1(2) + self.bin_num - self.preamble_bin)/self.zero_padding_ratio, 2^self.SF), ...
                  mod((pk_netid2(2) + self.bin_num - self.preamble_bin)/self.zero_padding_ratio, 2^self.SF)]
                ];
                
                % the goal is to extract payload_len from PHY header
                % header is in the first 8 symbols
                symbols = [];
                pk_list =[];
                if x > length(self.sig) - 8 * self.sample_num_per_symbol + 1
                    return;
                end
                for ii = 0:7
                    pk = self.dechirp(x + ii * self.sample_num_per_symbol);
                    pk_list = [pk_list; pk];
                    symbols = [symbols; mod((pk(2) + self.bin_num - self.preamble_bin)/self.zero_padding_ratio, 2^self.SF)];
                end
                if self.has_header
                    % get payload_len,crc_enable,cr
                    is_valid = self.check_header(symbols);
                    if ~is_valid
                        x = x + 7 * self.sample_num_per_symbol;
                        continue;
                    end
                end
                
                % number of symbols in the packet
                sym_num = self.calc_sym_num(self.Payload_len);
                
                % demodulate the rest LoRa data symbols
                if x > length(self.sig) - sym_num * self.sample_num_per_symbol + 1
                    return
                end
                for ii = 8:sym_num - 1
                    pk = self.dechirp(x + ii * self.sample_num_per_symbol);
                    pk_list = [pk_list; pk];
                    symbols = [symbols; mod((pk(2) + self.bin_num - self.preamble_bin)/self.zero_padding_ratio, 2^self.SF)];
                end  
                % compensate CFO drift
                symbols = self.dynamic_compensation(symbols);
                
                %                 x = x + sym_num * self.sample_num_per_symbol;
                
                symbols_m = [symbols_m mod(round(symbols), 2^self.SF)]; % symbols的值只能在[0,2^sf-1]之间
                cfo_m = [cfo_m self.CFO];
            end
            
            if isempty(symbols_m)
                warning('No preamble detected!');
            end
        end

        function x = pre_detect(self, start_idx)
            % Detect new form preamble 
            ii = start_idx;
            pk_bin_list = [];
            tmp_pattern_bin_list = [];
            self.pattern_bin_list = [];
            while ii < length(self.sig)-self.sample_num_per_symbol*self.preamble_len
                pk0 = self.pre_dechirp(ii);
                if size(pk_bin_list,2) < self.preamble_len
                    % 如果现在没有存够 6 个preamble chirp，则继续
                    % 否则追加到 pk_bin_list 中，并保持 pk_bin_list 的长度为8
                    pk_bin_list = [pk_bin_list pk0(:,2)];
                else
                    % 如果已经存够 6 个 preamble chirp,则检查是否是Preamble，若不是则追加
                    for i = 1:length(self.pattern)
                        row = self.pattern(i);
                        pk_bin = pk_bin_list(row,i);
                        self.pattern_bin_list = [self.pattern_bin_list pk_bin];
                    end
                    % 计算按照给定的 pattern 取出的bin_list，差值是否小于 2*padding

                    tmp_pattern_bin_list = self.pattern_bin_list - self.Ref_pk_bin_list;
                    tmp_pattern_bin_list = abs(tmp_pattern_bin_list);
                    found = false;
                    bin_diff = diff(tmp_pattern_bin_list);
                    if any(abs(bin_diff) > self.zero_padding_ratio*2) 
                        found = false;
                    else
                        found = true;
                    end
                    if found == true
                        warning('preamble detected!');
                        x = ii;
                        return;
                    else
                        pk_bin_list(:,1) = [];
                        pk_bin_list = [pk_bin_list pk0(:,2)];
                        self.pattern_bin_list = [];
                    end                    
                end
                ii = ii + self.sample_num_per_symbol;
            end
            x = -1;
        end

        function x_sync = pre_sync(self, x)
            % 功能：根据 pk_bin_list 求CFO和STO
            % inputs：
            %       pk_bin_list : {preamble + SFD} FFT bin
            
            up_pk_bin = 0;
            down_pk_bin = 0;
         
            if all(self.pattern == 1)
                up_pk_bin = mean(self.pattern_bin_list);
            elseif all(self.pattern == 2)
                down_pk_bin = mean(self.pattern_bin_list);
            else
                up_cnt = 0;
                down_cnt = 0;
                for i = 1:length(self.pattern)
                    if self.pattern(i) == 1
                        up_pk_bin = up_pk_bin + self.pattern_bin_list(i);
                        up_cnt = up_cnt + 1;
                    else
                        down_pk_bin = down_pk_bin + self.pattern_bin_list(i);
                        down_cnt = down_cnt + 1;
                    end
                end
                up_pk_bin = up_pk_bin/up_cnt;
                down_pk_bin = down_pk_bin/down_cnt;
            end
 
            cfo_offset = mod((up_pk_bin+down_pk_bin),self.fft_len*2)/2;
            
            % get -BW/2 corresponding fft bin
            fft_resolution = 1/(2*self.Ts)/self.zero_padding_ratio;
            Ref_down_bin = 3*self.BW/2/fft_resolution;
            
            % 根据 Bernier的论文，这种方式估计的cfo的范围在[−BW/4,BW/4]
%             if abs(cfo_offset) > self.bin_num/4
%                 if cfo_offset < 0
%                     cfo_offset = cfo_offset + self.bin_num/2;
%                 else
%                     cfo_offset = cfo_offset - self.bin_num/2;
%                 end
%             end
            
            self.preamble_bin = cfo_offset;

            sto_offset = round((down_pk_bin-cfo_offset-Ref_down_bin)/self.zero_padding_ratio);
            
%             if self.preamble_bin > self.bin_num / 2
%                 self.CFO = (self.preamble_bin-self.bin_num)*self.BW/self.bin_num;
%             else
%                 self.CFO = (self.preamble_bin)*self.BW/self.bin_num;
%             end
             
            % set x to the start of data symbols
%             pku = self.dechirp(x-2*self.sample_num_per_symbol+sto_offset);
%             pkd = self.dechirp(x-2*self.sample_num_per_symbol+sto_offset, false);
%             if abs(pku(1)) > abs(pkd(1))
%                 % current symbol is the first downchirp
%                 x_sync = x + round(2.25*self.sample_num_per_symbol) + sto_offset;
%             else
%                 % current symbol is the second downchirp
%                 x_sync = x + round(1.25*self.sample_num_per_symbol) + sto_offset;
%             end
            x_sync = x + round(2.25*self.sample_num_per_symbol) + sto_offset;
            if false
                sync = self.sig(x_sync:x_sync+10*self.sample_num_per_symbol);
                LoRaPHY.plot_timefrequency(sync,self.Fs,self.SF,self.BW);
                title("对齐之后，此时采样窗口的地方（位于第一个payload chirp)");
            end
        end

        function x = detect(self, start_idx)
            % Detect preamble
            %
            % input:
            %     start_idx: Start index for detection
            % output:
            %     x: Before index x, a preamble is detected.
            %        x = -1 if no preamble detected
            
            ii = start_idx;
            pk_bin_list = []; % 用来存放 preamble 的peak bin
            while ii < length(self.sig)-self.sample_num_per_symbol*self.preamble_len
                % search preamble_len-1 basic upchirps
                if length(pk_bin_list) == self.preamble_len - 1
                    % preamble detected
                    % coarse alignment: first shift the up peak to position 0
                    % current sampling frequency = 2 * bandwidth
            %                     x = ii - round((pk_bin_list(end)-1)/self.zero_padding_ratio*2); % 修改
                    x = ii;
                    return;
                end
                pk0 = self.dechirp(ii);
                if ~isempty(pk_bin_list)
                    bin_diff = mod(pk_bin_list(end)-pk0(2), self.bin_num);
                    if bin_diff > self.bin_num/2
                        bin_diff = self.bin_num - bin_diff; %考虑 bin_diff是负数的情况
                    end
                    if bin_diff <= self.zero_padding_ratio
                        pk_bin_list = [pk_bin_list; pk0(2)];
                    else
                        pk_bin_list = pk0(2);
                    end
                else
                    pk_bin_list = pk0(2);
                end
                ii = ii + self.sample_num_per_symbol;
            end
            x = -1;
        end
    
        
        function x_sync = sync(self, x)
            % sync  Packet synchronization
            %
            % input:
            %     x: Start index for synchronization
            % output:
            %     x_sync: Index after up-down alignment
            
            % find downchirp
            found = false;
            while x < length(self.sig) - self.sample_num_per_symbol
                up_peak = self.dechirp(x);
                down_peak = self.dechirp(x, false);
                if abs(down_peak(1)) > abs(up_peak(1))
                    % downchirp detected
                    found = true;
                end
%                 LoRaPHY.plot_timefrequency(self.sig(x:x+10*self.sample_num_per_symbol),self.Fs,self.SF,self.BW);
%                 title("sync定位到地方");
                x = x + self.sample_num_per_symbol;
                if found
                    break;
                end
            end
            
            if ~found
                return;
            end
            
            % Up-Down Alignment
            % NOTE preamble_len >= 6
            % NOTE there are two NETID chirps between preamble and SFD
            % NOTE `detect` has already shifted the up peak to position 0
            % NOTE bin-1 
            pkd = self.dechirp(x, false);
%             if pkd(2) > self.bin_num / 2
%                 to = round((pkd(2)-1-self.bin_num)/self.zero_padding_ratio);
%             else
%                 to = round((pkd(2)-1)/self.zero_padding_ratio);
%             end
%             x = x + to; 
            % 计算 SFD 的频率偏移: cfo + to;
            if pkd(2) > self.bin_num / 2
                sfd_offset = pkd(2)-self.bin_num-1;
            else 
                sfd_offset = pkd(2)-1;
            end
            
            % set preamble reference bin for CFO elimination
            % x - 5*self.sample_num_per_symbol 将指针指向preamble的base upchirp处
            pku = self.dechirp(x - 5*self.sample_num_per_symbol);
            
            % 计算 preamble 的频率偏移: cfo - to;
            if pku(2) > self.bin_num / 2
                preamble_offset = pku(2)-self.bin_num-1;
            else
                preamble_offset = pku(2)-1;
            end
            
            % 联立 sfd 和 preamble 的计算结果，得到 cfo
            cfo_offset = (sfd_offset+preamble_offset)/2;
            
            % 根据 Bernier的论文，这种方式估计的cfo的范围在[−BW/4,BW/4]
            if abs(cfo_offset) > self.bin_num/4
                if cfo_offset < 0
                    cfo_offset = cfo_offset + self.bin_num/2;
                else
                    cfo_offset = cfo_offset - self.bin_num/2;
                end
            end
            
            self.preamble_bin = cfo_offset;

            sto_offset = round((sfd_offset-cfo_offset)/self.zero_padding_ratio*2);
            
            if self.preamble_bin > self.bin_num / 2
                self.CFO = (self.preamble_bin-self.bin_num)*self.BW/self.bin_num;
            else
                self.CFO = (self.preamble_bin)*self.BW/self.bin_num;
            end
             
            % set x to the start of data symbols
            pku = self.dechirp(x-self.sample_num_per_symbol+sto_offset);
            pkd = self.dechirp(x-self.sample_num_per_symbol+sto_offset, false);
            if abs(pku(1)) > abs(pkd(1))
                % current symbol is the first downchirp
                x_sync = x + round(2.25*self.sample_num_per_symbol) + sto_offset;
            else
                % current symbol is the second downchirp
                x_sync = x + round(1.25*self.sample_num_per_symbol) + sto_offset;
            end
            if self.is_debug
                sync = self.sig(x_sync:x_sync+10*self.sample_num_per_symbol);
                LoRaPHY.plot_timefrequency(sync,self.Fs,self.SF,self.BW);
                title("对齐之后，此时采样窗口的地方（位于第一个payload chirp)");
            end
        end
        
        function is_valid = check_header(self, data)
            % check_header Parse LoRa PHY header and set parameters
            %
            % input:
            %     data: An eight elements vector containing header symbols
            % output:
            %     is_valid: `true` if the header is valid
            %               `false` if the header is invalid

            % compesate CFO drift
            symbols = self.dynamic_compensation(data);

            % gray coding
            symbols_g = self.gray_coding(symbols);

            %deinterleave
            codewords = self.diag_deinterleave(symbols_g(1:8), self.SF - 2);

            %解析 header
            nibbles = self.hamming_decode(codewords, 8);

            self.Payload_len = double(nibbles(1) * 16 + nibbles(2));
            self.CRC = double(bitand(nibbles(3), 1));
            self.CR = double(bitshift(nibbles(3), -1));

            % we only calculate header checksum on the first three nibbles
            % the valid header checksum is considered to be 5 bits
            % other 3 bits require further reverse engineering
            header_checksum = [bitand(nibbles(4), 1); de2bi(nibbles(5), 4, 'left-msb')'];
            header_checksum_calc = self.Header_checksum_matrix * gf(reshape(de2bi(nibbles(1:3), 4, 'left-msb')', [], 1));
            if any(header_checksum ~= header_checksum_calc)
                warning('Invalid header checksum!');
                is_valid = 0;
            else
                is_valid = 1;
            end
        
        end 
        %------------------------------------------------------------------------
        %
        % 以下是 decode 的代码
        %
        % ------------------------------------------------------------------------
        function [data_m, checksum_m] = decode(self, symbols_m)
            % decode  Decode data from symbols
            %
            % input:
            %     symbols_m: A matrix of symbols to be decoded. Each column
            %                vector represents demodulated symbols from a
            %                LoRa packet.
            % output:
            %     data_m: A matrix of bytes representing the decoding
            %             result of `symbols_m`. The last two bytes are the
            %             decoded CRC checksum if CRC is enabled.
            %     checksum_m: A vector of checksum based on the decoded
            %                 payload.Checksum is empty if CRC is disabled.

            data_m = [];
            checksum_m = [];

            for pkt_num = 1:size(symbols_m, 2)
                % gray coding
                symbols_g = self.gray_coding(symbols_m(:, pkt_num));

                % deinterleave
                codewords = self.diag_deinterleave(symbols_g(1:8), self.SF-2);
                if ~self.has_header
                    nibbles = self.hamming_decode(codewords, 8);
                else
                    % parse header
                    nibbles = self.hamming_decode(codewords, 8);
                    self.Payload_len = double(nibbles(1)*16 + nibbles(2));
                    self.CRC = double(bitand(nibbles(3), 1));
                    self.CR = double(bitshift(nibbles(3), -1));
                    % we only calculate header checksum on the first three nibbles
                    % the valid header checksum is considered to be 5 bits
                    % other 3 bits require further reverse engineering
                    header_checksum = [bitand(nibbles(4), 1); de2bi(nibbles(5), 4, 'left-msb')'];
                    header_checksum_calc = self.Header_checksum_matrix * gf(reshape(de2bi(nibbles(1:3), 4, 'left-msb')', [], 1));
                    if any(header_checksum ~= header_checksum_calc)
                        error('Invalid header checksum!');
                    end
                    nibbles = nibbles(6:end);
                end
                rdd = self.CR + 4;
                for ii = 9:rdd:length(symbols_g)-rdd+1
                    codewords = self.diag_deinterleave(symbols_g(ii:ii+rdd-1), self.SF-2*self.LDR);
                    % hamming decode
                    nibbles = [nibbles; self.hamming_decode(codewords, rdd)];
                end

                % combine nibbles to bytes
                bytes = uint8(zeros(min(255, floor(length(nibbles)/2)), 1));
                for ii = 1:length(bytes)
                    bytes(ii) = bitor(uint8(nibbles(2*ii-1)), 16*uint8(nibbles(2*ii)));
                end

                % dewhitening
                len = self.Payload_len;
                if self.CRC
                    % The last 2 bytes are CRC16 checkcum
                    data = [self.dewhiten(bytes(1:len)); bytes(len+1:len+2)];
                    % Calculate CRC checksum
                    checksum = self.calc_crc(data(1:len));
                else
                    data = self.dewhiten(bytes(1:len));
                    checksum = [];
                end
                data_m = [data_m data];
                checksum_m = [checksum_m checksum];
            end
        end

        function symbols = dynamic_compensation(self, data)
            % dynamic_compensation  Compensate bin drift
            %
            % input:
            %     data: Symbols with bin drift
            % output:
            %     symbols: Symbols after bin calibration

            % compensate the bin drift caused by Sampling Frequency Offset (SFO)
            % From Demodulation to Decoding equation(5)
            sfo_drift = (1 + (1:length(data))') * 2^self.SF * self.CFO / self.C_freq;
            
            symbols = mod(data - sfo_drift, 2^self.SF);

            if self.LDR
                bin_offset = 0;
                v_last = 1;
                
                for i = 1:length(symbols)
                    v = symbols(i);
                    bin_delta = mod(v - v_last, 4); %去掉后两个 bit 的数据
                    if bin_delta < 2
                        bin_offset = bin_offset - bin_delta;
                    else
                        bin_offset = bin_offset - bin_delta + 4;
                    end
                    v_last = v;
                    symbols(i) = mod(v + bin_offset, 2^self.SF);
                end
            end
        end

        function symbols = gray_coding(self, din)
            % gray_coding  Gray coding
            %              `gray_coding` is used in the DECODING process
            %
            % input:
            %     data: Symbols with bin drift
            % output:
            %     symbols: Symbols after bin calibration

            din(1:8) = floor(din(1:8)/4);
            if self.LDR
                din(9:end) = floor(din(9:end)/4);
            else
                din(9:end) = mod(din(9:end)-1, 2^self.SF);
            end
            s = uint16(din);
            symbols = bitxor(s, bitshift(s, -1));
            self.print_bin("Gray Coding", symbols, self.SF);
        end

        function codewords = diag_deinterleave(self, symbols, ppm)
            % diag_deinterleave  Diagonal deinterleaving
            %
            % input:
            %     symbols: Symbols after gray coding
            %     ppm: Size with parity bits
            % output:
            %     codewords: Codewords after deinterleaving

            b = de2bi(symbols, double(ppm), 'left-msb');
            codewords = flipud(bi2de(cell2mat(arrayfun(@(x) ...
                circshift(b(x,:), [1 1-x]), (1:length(symbols))', 'un', 0))'));
            self.print_bin("Deinterleave", codewords);
        end

        function nibbles = hamming_decode(self, codewords, rdd)
            % hamming_decode  Hamming Decoding
            %
            % input:
            %     codewords: Codewords after deinterleaving
            %     rdd: Bits with redundancy
            % output:
            %     nibbles: Nibbles after hamming decoding

            p1 = LoRaPHY.bit_reduce(@bitxor, codewords, [8 4 3 1]);
            p2 = LoRaPHY.bit_reduce(@bitxor, codewords, [7 4 2 1]);
            p3 = LoRaPHY.bit_reduce(@bitxor, codewords, [5 3 2 1]);
            p4 = LoRaPHY.bit_reduce(@bitxor, codewords, [5 4 3 2 1]);
            p5 = LoRaPHY.bit_reduce(@bitxor, codewords, [6 4 3 2]);
            function pf = parity_fix(p)
                switch p
                    case 3 % 011 wrong b3
                        pf = 4;
                    case 5 % 101 wrong b4
                        pf = 8;
                    case 6 % 110 wrong b1
                        pf = 1;
                    case 7 % 111 wrong b2
                        pf = 2;
                    otherwise
                        pf = 0;
                end
            end
            if self.hamming_decoding_en
                switch rdd
                    % TODO report parity error
                    case {5, 6}
                        nibbles = mod(codewords, 16);
                    case {7, 8}
                        parity = p2*4+p3*2+p5;
                        pf = arrayfun(@parity_fix, parity);
                        codewords = bitxor(codewords, uint16(pf));
                        nibbles = mod(codewords, 16);
                    otherwise
                        % THIS CASE SHOULD NOT HAPPEN
                        error('Invalid Code Rate!');
                end
            else
                nibbles = mod(codewords, 16);
            end
            self.print_bin("Hamming Decode", codewords);
        end

        function bytes_w = dewhiten(self, bytes)
            % dewhiten  Data Dewhitening
            %
            % input:
            %     bytes: Bytes after deinterleaving
            % output:
            %     bytes_w: Bytes after dewhitening

            len = length(bytes);
            bytes_w = bitxor(uint8(bytes(1:len)), self.Whitening_seq(1:len));
            self.print_bin("Dewhiten", bytes_w);
        end

        function pk = pre_dechirp(self, x)
            % 该函数用于 preamble 部分的dechirp，固定用两个解调窗口进行dechirp
            % dechirp 使用的 base up-/down-chirp 与一个解调窗口的相比，k不变，Ts*2

            % output : [bin_value  bin_index ]
            %          其中第一行使用 downchirp 进行解码
            %             第二行使用 upchirp 进行解码
            %              bin_value = -1 代表是没有peak
            upchirp = LoRaPHY.chirp(true,self.SF+2,2*self.BW,self.Fs,0);
            downchirp = LoRaPHY.chirp(false,self.SF+2,2*self.BW,self.Fs,0);
           
            ft_up = fft(self.sig(x:x+self.sample_num_per_symbol*2-1).*upchirp, self.fft_len*2);
            ft_up_ = abs(ft_up);
            ft_down = fft(self.sig(x:x+self.sample_num_per_symbol*2-1).*downchirp, self.fft_len*2);
            ft_down_ = abs(ft_down);
            % --- DEBUG ---
            if false
                signal = self.sig(x:x+self.sample_num_per_symbol*2-1);
                LoRaPHY.plot_timefrequency(signal,self.Fs,self.SF,self.BW);
                title("解调窗口的信号时频图");
%                 LoRaPHY.plot_timefrequency(upchirp,self.Fs,self.SF,self.BW);
%                 title("dechirp使用的base chirp的时频图");
%                 LoRaPHY.plot_timefrequency(downchirp,self.Fs,self.SF,self.BW);
%                 title("dechirp使用的base chirp的时频图");
                
                figure;
                subplot(211);
                plot(ft_up_);
                title("preamble使用两个解调窗口的FFT结果图(使用upchirp)");
                xlabel('bin');
                subplot(212);
                plot(ft_down_);
                title("preamble使用两个解调窗口的FFT结果图(使用downchirp)");
                xlabel('bin');
            end
            is_Outlier_up = LoRaPHY.hasOutlier([ft_up_ (1:self.fft_len*2).']);
            is_Outlier_down = LoRaPHY.hasOutlier([ft_down_ (1:self.fft_len*2).']);
            if is_Outlier_up == is_Outlier_down
                % 如果在down/up-demodulation window 中均存在离群点或均不存在离群点时，
                % 则此时存在解调窗口存在 upchirp和downchirp，那么则各自取peak
                pk_up = LoRaPHY.topn([ft_up_ (1:self.fft_len*2).'], 1);
                pk_down = LoRaPHY.topn([ft_down_ (1:self.fft_len*2).'], 1);
            end
            if is_Outlier_up && ~is_Outlier_down
                % 如果在up-demodulation window 中存在离群点
                % 在 down-demodulation window 中不存在离群点
                % 则 此时解调窗口存在两个downchirp，
                % 如果则合法的peak应该出现在[3B/2,2B]处
                pk_down = [-1 -1]; 
                pk_up = LoRaPHY.topn([ft_up_(3*self.fft_len/2:self.fft_len*2) (3*self.fft_len/2:self.fft_len*2).'], 1);
            end
            if ~is_Outlier_up && is_Outlier_down
               % 如果在down-demodulation window 中存在离群点
               % 在 up-demodulation window 中不存在离群点
               % 则 此时解调窗口存在两个upchirp，
               % 如果则合法的peak应该出现在[1,B/2]处
               pk_up = [-1 -1]; 
               pk_down = LoRaPHY.topn([ft_down_(1:self.fft_len/2) (1:self.fft_len/2).'], 1);
            end
            pk = [pk_down;pk_up];
        end

        function pk = dechirp(self, x, is_up)
            % dechirp  Apply dechirping on the symbol starts from index x
            %
            % input:
            %     x: Start index of a symbol
            %     is_up: `true` if applying up-chirp dechirping
            %            `false` if applying down-chirp dechirping
            % output:
            %     pk: Peak in FFT results of dechirping
            %     pk = (height, index)
        
            if nargin == 3 && ~is_up
                c = self.Ref_upchirp;
            else
                c = self.Ref_downchirp;
            end
            ft = fft(self.sig(x:x+self.sample_num_per_symbol-1).*c, self.fft_len);
            ft_ = abs(ft(1:self.bin_num)) + abs(ft(self.fft_len-self.bin_num+1:self.fft_len)); %合并两个peak,Coarse-grained Phase Alignment (CPA) 
%             if self.is_debug
%                 fft_bin = (0:self.bin_num-1);
%                 figure;
%                 title("symbol的fft结果图")
%                 plot(fft_bin,ft_);
%                 xlabel('bin')
%             end
            pk = LoRaPHY.topn([ft_ (1:self.bin_num).'], 1);
        end
        
        function print_bin(self, flag, vec, size)
            if self.is_debug
                if nargin == 3
                    size = 8;
                end
                len = length(vec);
                fprintf("%s:\n", flag);
                for i = 1:len
                    fprintf("%s\n", dec2bin(round(vec(i)), size));
                end
                fprintf("\n");
            end
        end
     
    end

    methods(Static)

        function b = bit_reduce(fn, w, pos)
            b = bitget(w, pos(1));
            for i = 2:length(pos)
                b = fn(b, bitget(w, pos(i)));
            end
        end

        function w = word_reduce(fn, ws)
            w = ws(1);
            for i = 2:length(ws)
                w = fn(w, ws(i));
            end
        end

        function print_payload(symbols)
            % print payload using symbols array
            symbols = symbols(1:end-2,:);
            ss = char(symbols);
            columnStrings = arrayfun(@(col) ss(:, col)', 1:size(ss, 2), 'UniformOutput', false);
            for i = 1:numel(columnStrings)
                disp(columnStrings{i});
            end
        end

        function y = chirp(is_up,sf,bw,fs,h,cfo,tdelta,tscale)
            % chirp  Generate a LoRa chirp symbol
            %
            % input:
            %     is_up: `true` if constructing an up-chirp
            %            `false` if constructing a down-chirp
            %     sf: Spreading Factor
            %     bw: Bandwidth
            %     fs: Sampling Frequency
            %     h: Start frequency offset (0 to 2^SF-1)
            %     cfo: Carrier Frequency Offset
            %     tdelta: Time offset (0 to 1/fs)
            %     tscale: Scaling the sampling frequency
            % output:
            %     y: Generated LoRa symbol
            if nargin < 8
                tscale = 1;
            end
            if nargin < 7
                tdelta = 0;
            end
            if nargin < 6
                cfo = 0;
            end
            N = 2^sf;
            T = N/bw;
            sample_per_symbol = round(fs * T);
            h_orig = h;
            h = round(h);
            cfo = cfo + (h_orig - h) * bw / N;
            
            if is_up
                k = bw/T;
                f0 = -bw/2 + cfo;
            else
                k = -bw/T;
                f0 = bw/2 + cfo;
            end
            t1 = (0:sample_per_symbol*(N-h)/N)/fs*tscale + tdelta; % first segment chirp
            snum = length(t1);
            c1 = exp(1j*2*pi*(t1.*(f0 + bw*h/N + k/2*t1)));
            if snum == 0
                phi = 0;
            else 
                phi = angle(c1(snum));
            end
            t2 = (0:sample_per_symbol*h/N-1)/fs + tdelta; % second segment chirp
            c2 = exp(1j*(phi + 2*pi*(t2.*(f0 + k/2*t2))));
            
            y = cat(2,c1(1:snum-1),c2).';
        end
        
        function plot_timefrequency(signal,fs,sf,bw)
            % 画出的时频图
            win_length = 2^(sf-2);
            Ts = 2^(sf)/bw;
            nfft = fs * Ts;
            s = spectrogram(signal,win_length,round(win_length*0.8),nfft);
            s = fftshift(s,1);
            figure;
            x = 0:1/fs:2^sf/bw;
            y = -bw/2:bw/2;
            imagesc(x,y,abs(s));
            title('Spectrogram');
            xlabel('Time/(ms)');
            ylabel('Frequency/(Hz)');
            set(gca,'ydir','normal');
%             [s,f,t] = spectrogram(data,256,128,nfft,fs);
%             figure;
%             waterplot(s,f,t);
%             title('Spectrogram');
%             xlabel('Time');
%             ylabel('Frequency');
%             set(gca,'ydir','normal');
        
        end

        function waterplot(s,f,t)
            % Waterfall plot of spectrogram
            waterfall(f,t,abs(s)'.^2)
            set(gca,XDir="reverse",View=[30 50])
            xlabel("Frequency (Hz)")
            ylabel("Time (s)")
        end

        function iq_values = read_file(filename)
            %读入文件
            f = fopen(filename,'rb');
            values = fread(f,[2, Inf],'float');
            fclose(f);
            iq_values = values(1,:) + values(2,:)*1i;
            [r, c] = size (iq_values);
            iq_values = reshape (iq_values, c, r);
        end

        function v = write(data, filename)
        
            % usage: write(data, filename)
            %
            %  open filename and write data to it
            %  Format is interleaved float IQ e.g. each
            %  I,Q 32-bit float IQIQIQ....
            %  This is compatible with read_complex_binary()
            %
            
            m = nargchk (2,2,nargin);
            if (m)
            usage (m);
            end
            
            f = fopen (filename, 'wb');
            if (f < 0)
            v = 0;
            else
            re = real(data);
            im = imag(data);
            re = re(:)';
            im = im(:)';
            y = [re;im];
            y = y(:);
            v = fwrite (f, y, 'float');
            fclose (f);
            end
        end
        function outlier = hasOutlier(pks)
            
            % 计算均值和方差
            data = abs(pks(:,1));
            mean_val = mean(data);
            var_val = std(data);
            threshold = mean_val + 5*var_val;
            outlier_idx = find(data > threshold); % 返回离群点的idx
            if isempty(outlier_idx)
                outlier = false;
            else
                outlier = true;
            end
                

        end

        function y = topn(pks, n, padding, th)
        %             input:
        %                 fft_res : fft的结果
        %                 n : 返回前 n 个振幅最大的峰值
        %                 padding：一个可选的参数，用于指定是否需要对结果进行填充以确保其长度为 n;默认值为 false
        %                 th：一个可选的阈值，用于指定需要在振幅小于该值的位置处停止选择峰值。默认情况下，函数会选择所有振幅最大的峰值。
        %             output:
        %                 函数的输出是一个 n×2 的数组，其中每一行都包含一个峰值的振幅和对应的index。
        %                 如果 padding 参数为 true，则输出数组的长度可能小于 n。如果提供了 th 参数，则输出数组中可能不包含所有振幅大于阈值的峰值
            [y, p] = sort(abs(pks(:,1)), 'descend');
            if nargin == 1
                return;
            end
            nn = min(n, size(pks, 1));
            if nargin >= 3 && padding
                y = [pks(p(1:nn), :); zeros(n-nn, size(pks, 2))];
            else
                y = pks(p(1:nn), :);
            end
            if nargin == 4
                ii = 1;
                while ii <= size(y,1)
                    if abs(y(ii,1)) < th
                        break;
                    end
                    ii = ii + 1;
                end
                y = y(1:ii-1, :);
            end
        end

        function peak_index = fft_by_symbol(signal,sf,os_factor,is_plot)
        %计算一个symbol的fft结果
        % input:
        %                   is_plot : 是否要画出 fft结果图
        %             output:
        %                   peak_index: Peak in FFT results of dechirping
        %                   peak_index = (height, index)
            nfft = length(signal);
            fft_res = fft(signal);
            fft_res_ = abs(fft_res);
            if nargin >=3
                peak_index = LoRaPHY.topn([fft_res_ (1:nfft).'],1);
            end
            if nargin == 4 && is_plot
                fft_bin = (0:nfft-1);
                figure;
                title("symbol的fft结果图")
                plot(fft_bin,fft_res_);
                xlabel('bin')
            end
        end
        
        function varargout = dechirp_(sig,fs,bw,sf,is_up,dw_n,is_debug,zero_padding)
%             input:
%                 dw_n: the number of demodulation window
%                 is_up: `true` if constructing an up-chirp
%                        `false` if constructing a down-chirp
%                 zero_padding: the scale of fft 

            if nargin >= 5 && ~is_up
                c = LoRaPHY.chirp(true,sf,bw,fs,0,0);
            else
                c = LoRaPHY.chirp(false,sf,bw,fs,0,0);
            end
            if nargin >= 6 && dw_n ~= 1
                c = repmat(c,dw_n,1);
            else 
                dw_n = 1;
            end
            if nargin < 7
                is_debug = false;
            end
            if nargin < 8
                zero_padding = 1;
            end
            N = 2^(sf);
            Ts = N/bw;
            fft_len = fs * Ts * dw_n * zero_padding;
            bin_num = N * dw_n * zero_padding;
%             ft = fft(sig,fft_len);
            ft = fft(sig.*c,fft_len);
%             ft_ = abs(ft(1:bin_num));
            ft_ = abs(ft(1:bin_num)) + abs(ft(fft_len-bin_num+1:fft_len));
            if is_debug
                x_lable = (0:bin_num-1);
                figure;
                title("symbol的fft结果图")
                plot(x_lable,ft_);
                xlabel('bin')
            end
            pk = LoRaPHY.topn([ft_ (1:bin_num).'], 1);

            varargout{1} = pk;
            if nargout == 2
                varargout{2} = ft;
            elseif nargout == 3
                varargout{2} = ft;
                varargout{3} = ft_;
            end
        end

        function plot_top_n(A,n)
        % 该函数将数组A画出，并标记处前n个最大的值
        % 注意，A应该是实数
            [sortedValues, sortedIndices] = sort(A, 'descend');
            topNIndices = sortedIndices(1:n);
            figure;
            plot(A);
            hold on;
            scatter(topNIndices, A(topNIndices), 'r', 'filled');
        end
    end    
end