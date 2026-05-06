classdef OFDM
    methods(Static)
        function S = generateOFDMframe(params)
            Nf = params.numFrames;
            N = params.numCarriers;
            blocks_per_frame = params.preamblePeriodicity;
            data_blocks_per_frame =  blocks_per_frame - 1;

            if params.pilots
                pilot_indices = OFDM.pilot_indices(params);
            else
                pilot_indices = [];
            end
            
            Ns = Nf*blocks_per_frame;
            S = zeros([N Ns]);
            
            guard_indices = [params.guardTones, params.DCSubcarrier] + N/2 ;
            
            excluded = unique([guard_indices(:); pilot_indices(:)]);
            data_indices = setdiff((1:N).', excluded, 'stable');
            
            % Place preambles at the start of each frame (column 1 of each frame block)
            preamble = OFDM.generatePreamble(params);
            preamble_col_indices = 1 : blocks_per_frame : Ns;
            S(:, preamble_col_indices) = repmat(preamble, 1, numel(preamble_col_indices));
            
            data_col_indices = setdiff(1:Ns, preamble_col_indices, 'stable');

            % Fill by columns
            Sdata = OFDM.generateSymbols(params, length(data_indices)* Nf*data_blocks_per_frame);
            S(data_indices, data_col_indices) = reshape(Sdata, numel(data_indices),[]);

            if params.pilots
                S = OFDM.add_pilots(S, params);
            end
        end

        function data_frame_indices = data_frame_indices(params)
            Nf = params.numFrames;
            blocks_per_frame = params.preamblePeriodicity;
            Ns = Nf*blocks_per_frame;
            preamble_col_indices = 1 : blocks_per_frame : Ns;
            data_frame_indices = setdiff(1:Ns, preamble_col_indices, 'stable');
        end

        function R = demodulate(r, params)
            N = params.numCarriers;
            Lcp = params.prefixLength;
            r = reshape(r, N+Lcp, []);
            r = OFDM.remove_prefix(r, params);  % Remove CP before FFT
            R = OFDM.DFT(r, params)/sqrt(N);
        end

        function h = channelEstimate(R, params)
            pilot_indices = OFDM.pilot_indices(params);
            data_frame_indices = OFDM.data_frame_indices(params);
            N = params.numCarriers;
            numSymbols = length(data_frame_indices);
            ref = OFDM.getPilots(params);
            Rp = R(pilot_indices, data_frame_indices);  % Already correctly scaled from demodulate
            
            if lower(params.estimationType) == "zf"
                H = zeros(N, numSymbols);
                H(pilot_indices, :) = bsxfun(@rdivide, Rp, ref);
                h = OFDM.IDFT(H, params)*sqrt(N);
            elseif lower(params.estimationType) == "zf_linear_interpolate"
                % First, get channel at pilot locations
                H_pilots = bsxfun(@rdivide, Rp, ref);  % num_pilots x numSymbols
                
                % For each symbol, linearly interpolate H to all subcarriers
                H = zeros(N, numSymbols);
                for sym = 1:numSymbols
                    % Interpolate from pilot locations to all subcarriers
                    H(:, sym) = interp1(pilot_indices, H_pilots(:, sym), (1:N).', 'linear', 'extrap');
                end
                
                h = OFDM.IDFT(H, params);%*sqrt(N);

            elseif lower(params.estimationType) == "zf_fft_interpolate"
                % FFT-based interpolation: zero-pad in frequency domain
                H_pilots = bsxfun(@rdivide, Rp, ref);  % num_pilots x numSymbols
                
                % For each symbol, use FFT-based interpolation
                H = zeros(N, numSymbols);
                for sym = 1:numSymbols
                    % Create sparse frequency domain representation at pilot locations only
                    H_sparse = zeros(N, 1);
                    H_sparse(pilot_indices) = H_pilots(:, sym);
                    
                    % Transform to time domain
                    h_time = OFDM.IDFT(H_sparse, params);
                    
                    % FFT back to frequency domain (gives interpolated response at all subcarriers)
                    H(:, sym) = OFDM.DFT(h_time, params);
                end
                
                h = OFDM.IDFT(H, params)*sqrt(N);
            end
                
                h = OFDM.IDFT(H, params);
            end
        end
        
        function mse = channelMSE(h, hest, params)
            mse = 0;
            h = h.';
            N = params.numCarriers;
            h = [h ;zeros([N-length(h) 1])];
            for i = 1:size(hest, 2)
                mse = mse + mean(abs(h-hest(:, i)).^2);
            end

        end

        function s = modulate(S, params)
            N = params.numCarriers;
            s = OFDM.IDFT(S, params)*sqrt(N);
            s = OFDM.add_prefix(s, params);
            s = s(:);
        end
        function pilot_indices = pilot_indices(params) 
            N = params.numCarriers;         
            start_tone = min(params.upperGuardTones)-1;
            end_tone = max(params.lowerGuardTones)+1;
            pilot_indices = (start_tone:-params.pilotPeriodicity:end_tone) + N/2;
        end
        function Sp = add_pilots(S, params)
            idx = OFDM.pilot_indices(params);
            Sp = S;
            pilots = OFDM.getPilots(params);

            data_indices = OFDM.data_frame_indices(params);
            % pilots is already a vector (length = num_pilots)
            % Assign to all columns
            Sp(idx, data_indices) = repmat(pilots, 1, length(data_indices));
        end

        function pilots = getPilots(params)
           idx = OFDM.pilot_indices(params);
           pilot_val = ((1+1j)/sqrt(2));
           pilots = repmat(pilot_val, length(idx),1);
        end

        function scp = add_prefix(s, params)

            N  = params.numCarriers;
            Lcp = params.prefixLength;

            n = (-Lcp:-1).';
            idx = mod(N+n, N)+1;
            s_prefix = s(idx, :);
            scp = [s_prefix; s];
        end

        function s = remove_prefix(scp, params)
        %REMOVE_PREFIX Remove prefix from OFDM signal
        %
        %   s = ofdm.remove_prefix(scp, params)
        %
        %   INPUTS:
        %       scp    : Signal with prefix
        %       params : Structure with field prefixLength
        %
        %   OUTPUT:
        %       s : Prefix-removed signal
            Lcp = params.prefixLength;
            s = scp(Lcp+1:end,:);
        end

        function s = IDFT(F, params)
            % IDFT using IFFT with OFDM-style shift (DC centered input).
            % Input F is N x numSamples, ordered as [-N/2 ... N/2-1] along rows.
            % Output s is N x numSamples (time-domain), column-wise signals.
            N = params.numCarriers;
            s = ifft(ifftshift(F, 1), N, 1);
        end

        function F = DFT(s, params)
            % DFT using FFT with OFDM-style shift (DC centered output).
            % Input s is N x numSamples along rows/time.
            % Output F is N x numSamples, ordered as [-N/2 ... N/2-1] along rows.
            N = params.numCarriers;
            F = fftshift(fft(s, N, 1), 1);
        end


        %%%%% Schmidl Cox Preamble %%%%%
        function preamble = generatePreamble(params)
            N = params.numCarriers;

            preamble_indices = (1:2:N).';
            preamble = zeros(N, 1);

            % CAZAC Preamble
            % num_preamble_symbols = length(preamble_indices);
            % cazac_seq = OFDM.generateCAZAC(num_preamble_symbols); 
            % preamble(preamble_indices) = cazac_seq;

            % % QPSK preamble
            qpskParams = params;
            qpskParams.modulationType  = 'QAM';
            qpskParams.modulationOrder = 4;
            preamble(preamble_indices) = OFDM.generateSymbols(qpskParams, numel(preamble_indices));


            guard_indices = [params.guardTones, params.DCSubcarrier] + N/2 ;
            preamble(guard_indices) = 0;            
        end

        function cazac_seq = generateCAZAC(N)
            % Generate Zadoff-Chu (CAZAC) sequence of length N
            % u is the root (typically coprime to N)
            % cazac_seq(n) = exp(j * pi * u * n * (n+1) / N) for n = 0:N-1
            
            if N == 1
                cazac_seq = 1;
                return;
            end
            
            % Find a root u coprime to N
            u = 1;
            for candidate = 1:N
                if gcd(candidate, N) == 1
                    u = candidate;
                    break;
                end
            end
            
            n = (0:N-1).';
            cazac_seq = exp(1j * pi * u * n .* (n + 1) / N);
        end

        function symbols = generateSymbols(params, num_data_symbols)
            % Prefer correct spelling first
            M = params.modulationOrder;
            idx = randi([0, M-1], num_data_symbols, 1);

            switch upper(params.modulationType)
                case 'QAM'
                    symbols = qammod(idx, M, 'UnitAveragePower', true);

                case 'PSK'
                    % QPSK diagonal constellation when M=4:
                    symbols = pskmod(idx, M, pi/4);

            end
        end
        

        function [r, h] = channel(s, params)
            h = OFDM.createChannelRealization(params);
            r_full = conv(s, h);
            r = r_full(1:length(s));

            % r = s;
            r = awgn(r, params.SNR, 'measured');
            cfo = params.CFO;
            Ts = 1/params.sampleRate;
            Ns = length(r);
            n = 0:Ns-1;
            r = r.';
            r = r.*exp(1j*2*pi*cfo*n*Ts);
        end


        function h = createChannelRealization(params)
            tau = params.delays;
            pg = 10.^(params.pathGains/10);
            tau_max = max(tau);

            h = zeros([1, tau_max+1]);

            for i = 1: length(tau)
                h(tau(i)+1) = sqrt(pg(i))*(randn(1, 1)+1j*randn(1,1))/sqrt(2);
            end
        end

        function f_est = estimateCFO(r, params)
            N = params.numCarriers;
            Lcp = params.prefixLength;
            Np = params.preamblePeriodicity;
            Ts = 1/params.sampleRate;

            r = reshape(r, N + Lcp, []);
            rp = r(:,1:Np:end);

            f_est = zeros(1,size(rp, 2));
            for n = 1:size(rp,2)
                p =  0;
                for i = Lcp+1:Lcp+N/2
                    p = p + rp(i, n)*rp(i+N/2, n)';
                end
                phase = angle(p');
                f_est(n) = phase/(N*Ts*pi);
            end
        end

        function z = SCcorr(r, params)
            N = params.numCarriers;
            Lcp = params.prefixLength;
            Np = params.preamblePeriodicity;
            % r = reshape(r, N + Lcp, []);
            
            
            % rp = r(:,1:Np:end);
            % numPreambles = size(rp, 2);
            % win = N/2;
            % z = zeros(Lcp+1, numPreambles);
            % for col = 1:numPreambles
                
            %     x = [rp(:, col).'];
            %     for d = 1:length(x)-N;
            %         a = x(d : d + win- 1);
            %         b = x(d + N/2 : d + N/2+win - 1);
            %         P = sum(a .* conj(b));
            %         R = sum(abs(b).^2);
            %         if R > 0
            %             z(d, col) = abs(P)^2 / R^2;
            %         end
            %     end
            % end

            rp = r(:);
            z = zeros([1 length(rp)-N]);
            for i=1:length(rp)-N
                a = rp(i : i + N/2 -1);
                b = rp(i + N/2: i + N-1);
                P = sum( a.* conj(b));
                R = 0.5* (sum(abs(b).^2) + sum(abs(a).^2));
                z(i) = abs(P^2)/R^2;
                % z(i) = abs(P^2);
            end
        end

    end
end
