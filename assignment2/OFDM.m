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

        function data_indices = data_subcarrier_indices(params)
            N = params.numCarriers;
            if isfield(params, 'pilots') && params.pilots
                pilot_indices = OFDM.pilot_indices(params);
            else
                pilot_indices = [];
            end
            guard_indices = [params.guardTones, params.DCSubcarrier] + N/2 ;
            
            excluded = unique([guard_indices(:); pilot_indices(:)]);
            data_indices = setdiff((1:N).', excluded, 'stable');
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
            R = OFDM.DFT(r, params);%/sqrt(N);
        end

        function h = channelEstimate(R, params)
            pilot_indices = OFDM.pilot_indices(params);
            data_frame_indices = OFDM.data_frame_indices(params);
            N = params.numCarriers;
            Lcp = params.prefixLength;
            numSymbols = length(data_frame_indices);
            ref = OFDM.getPilots(params);
            Rp = R(pilot_indices, data_frame_indices);  % Already correctly scaled from demodulate
            
            estimationType = lower(string(params.estimationType));

            if estimationType == "zf"
                H = zeros(N, numSymbols);
                H(pilot_indices, :) = bsxfun(@rdivide, Rp, ref);
                h = OFDM.IDFT(H, params);

            elseif estimationType == "zf_linear_interpolate"
                % First, get channel at pilot locations
                H_pilots = bsxfun(@rdivide, Rp, ref);  % num_pilots x numSymbols
                
                % For each symbol, linearly interpolate H to all subcarriers
                H = zeros(N, numSymbols);
                for sym = 1:numSymbols
                    % Interpolate from pilot locations to all subcarriers
                    H(:, sym) = interp1(pilot_indices, H_pilots(:, sym), (1:N).', 'linear', 'extrap');
                end
                
                h = OFDM.IDFT(H, params);

            elseif lower(params.estimationType) == "zf_fft_interpolate" 
                H_pilots = bsxfun(@rdivide, Rp, ref);  % num_pilots x numSymbols
                H = zeros(N, numSymbols);
                [pilot_idx_sorted, sort_order] = sort(pilot_indices(:));
                query_idx = (1:N).';
                fft_grid_len = length(pilot_idx_sorted) * params.pilotPeriodicity;
                fft_grid = pilot_idx_sorted(1) + (0:fft_grid_len-1).';
                for sym = 1:numSymbols
                    H_fft_interp = interpft(H_pilots(sort_order, sym), fft_grid_len);
                    H(:, sym) = interp1(fft_grid, H_fft_interp, query_idx, 'linear', 'extrap');
                end
                h = OFDM.IDFT(H, params);
            elseif lower(params.estimationType) == "mls"
                h = zeros(Lcp, numSymbols);       
                for i = 1:numSymbols
                    h(:, i) = OFDM.mLSestimation(Rp(:, i), ref, params);
                end
                h = [h ; zeros([N-Lcp numSymbols])]/sqrt(N);
            
            elseif lower(params.estimationType) == "mls_informed"
                Nt = length(params.delays);

                h = zeros(N, numSymbols);       
                for i = 1:numSymbols
                    h(params.delays+1, i) = OFDM.mLSestimation_informed(Rp(:, i), ref, params);
                end
                h = h/sqrt(N);
            end
        end
        
        function h = mLSestimation(Yp, xp, params)
            Lcp = params.prefixLength;
            pilot_indices = OFDM.pilot_indices(params);
            F = OFDM.getDFTMatrix(params);
            Fp = F(pilot_indices, 1:Lcp);
            Xp = diag(xp);
            A = Xp*Fp;

            h = (A' * A ) \  A' * Yp;
        end
        function h = mLSestimation_informed(Yp, xp, params)
            Lcp = params.prefixLength;
            pilot_indices = OFDM.pilot_indices(params);
            F = OFDM.getDFTMatrix(params);
            Fp = F(pilot_indices, params.delays+1);
            Xp = diag(xp);
            A = Xp*Fp;

            h = (A' * A ) \  A' * Yp;
        end

        function mse = channelMSE(h, hest, params)
            mse = 0;
            N = params.numCarriers;
            estimationType = lower(string(params.estimationType));

            excluded_indices = [params.guardTones(:); params.DCSubcarrier(:)] + N/2;
            active_indices = setdiff((1:N).', excluded_indices, 'stable');

            H = OFDM.DFT(h(:), params);
            for i = 1:size(hest, 2)
                Hest = OFDM.DFT(hest(:, i), params);
                if estimationType == "zf"
                    pilot_indices = OFDM.pilot_indices(params);
                    Hcmp = H(pilot_indices);
                    Hestcmp = Hest(pilot_indices);
                else
                    Hcmp = H(active_indices);
                    Hestcmp = Hest(active_indices);
                end

                mse = mse + mean(abs(Hcmp - Hestcmp).^2);
            end
            mse = mse / size(hest,2);
        end

        function s = modulate(S, params)
            N = params.numCarriers;
            s = OFDM.IDFT(S, params);%*sqrt(N);
            s = OFDM.add_prefix(s, params);
            s = s(:);
        end
        function pilot_indices = pilot_indices(params) 
            N = params.numCarriers;
            if isfield(params, 'upperGuardTones') && isfield(params, 'lowerGuardTones')
                upper_guard_tones = params.upperGuardTones;
                lower_guard_tones = params.lowerGuardTones;
            else
                upper_guard_tones = params.guardTones(params.guardTones > 0);
                lower_guard_tones = params.guardTones(params.guardTones < 0);
            end

            start_tone = min(upper_guard_tones)-1;
            end_tone = max(lower_guard_tones)+1;
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

        function W = getDFTMatrix(params)
            % Generate DFT matrix for OFDM demodulation with DC-centered ordering
            % 
            % Returns: W - DFT matrix (N x N) with fftshift applied
            %   Usage: F = W * s, where s is time-domain vector
            %   Output F is DC-centered: [-N/2 ... N/2-1]
            
            N = params.numCarriers;
            
            % Standard DFT matrix (unshifted)
            n = (0:N-1)';
            k = (0:N-1);
            W_standard = exp(-1j * 2 * pi * n * k / N) / sqrt(N);
            
            % Apply fftshift to get DC-centered ordering
            % fftshift moves [0...N-1] to [-N/2...N/2-1]
            W = fftshift(W_standard, 1);
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
            h = h.';
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

        function plotChannelEstimationPerformance(snrs, mses, estimationTypes)
            % Plot channel estimation MSE vs SNR for different methods
            %
            % Usage:
            %   OFDM.plotChannelEstimationPerformance(snrs, mses, estimationTypes)
            %
            % Inputs:
            %   snrs: SNR values in dB (1 x numSNRs)
            %   mses: MSE matrix (numMethods x numSNRs)
            %   estimationTypes: String array of method names
            
            mses_db = 10*log10(mses);
            
            figure('Position', [100, 100, 1000, 600]);
            hold on;
            grid on;
            
            % Define colors and markers for different methods
            colors = lines(length(estimationTypes));
            markers = {'o', 's', '^', 'd', 'v'};
            
            % Plot each estimation method
            for i = 1:length(estimationTypes)
                marker_idx = mod(i-1, length(markers)) + 1;
                plot(snrs, mses_db(i, :), ...
                    'Color', colors(i, :), ...
                    'Marker', markers{marker_idx}, ...
                    'LineWidth', 2, ...
                    'MarkerSize', 8, ...
                    'DisplayName', estimationTypes(i));
            end
            
            % Labels and formatting
            xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('Channel Estimation MSE (dB)', 'FontSize', 12, 'FontWeight', 'bold');
            title('Channel Estimation Performance Comparison', 'FontSize', 14, 'FontWeight', 'bold');
            
            legend('Location', 'best', 'FontSize', 11, 'Interpreter', 'none');
            grid on;
            grid minor;
            
            % Set axis properties
            xlim([snrs(1), snrs(end)]);
            set(gca, 'FontSize', 11);
            
            hold off;
        end

    end
end
