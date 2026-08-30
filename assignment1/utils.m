classdef utils
methods(Static)
    function h = levinsonFilter(fd, Ts, p, N)
        % --- Theoretical autocorrelation ---
        m = 0:p;
        R = besselj(0, 2*pi*fd*m*Ts);   % Jakes autocorrelation
        R(1) = R(1) + 1e-6;
        
        % --- Levinson-Durbin ---
        [a, E] = levinson(R, p);  % a = [1 a1 ... ap]
        
        % --- Generate white complex Gaussian ---
        w = (randn(1, N) + 1j*randn(1, N))/sqrt(2);
        
        % --- Generate colored process ---
        h = filter(sqrt(E), a, w);
        
    end

    function [h, t] = smithsFilter(fd, Ts, N)
        df = 1/(N*Ts);
        f  = (-N/2:N/2-1) * df;

        % --- Jakes PSD ---
        S = zeros(1, N);
        idx = abs(f) <= fd;
        
        S(idx) = 1 ./ sqrt(1 - (f(idx)/fd).^2 + eps);  % avoid infinity
        
        % --- Generate complex white Gaussian noise ---
        W = (randn(1,N) + 1j*randn(1,N))/sqrt(2);
        
        % --- Shape spectrum ---
        H = sqrt(S) .* W;
        
        % --- Time domain fading process ---
        h = ifft(ifftshift(H));
        
        % --- Normalize power ---
        h = h / sqrt(mean(abs(h).^2));

        % --- Time axis ---
        t = (0:N-1)*Ts;
    end

    function [h, t] = SOS(fd, Ts, N0, N)
        t = (0:N-1)*Ts;
        
        % angles
        n = 1:N0;
        theta = pi*(n - 0.5)/N0;
        
        % Doppler frequencies
        fn = fd * cos(theta);
        
        % random phases
        phi = 2*pi*rand(1,N0);
        psi = 2*pi*rand(1,N0);
        
        % generate fading
        h = zeros(1,N);
        
        for k = 1:N0
            h = h + ...
                cos(2*pi*fn(k)*t + phi(k)) + ...
                1j*sin(2*pi*fn(k)*t + psi(k));
        end
        
        % normalize
        h = sqrt(2/N0) * h;
    end
    function [h, t] = channelRealization(pg, tau, Ts, Nfft)

        pg = pg / sum(pg); % normalize total power = 1
        
        delay_idx = round(tau / Ts);
        
        h = zeros(1, Nfft); % directly create padded version
        
        for i = 1:length(pg)
            w = (randn + 1j*randn)/sqrt(2); % CN(0,1)
            h(delay_idx(i) + 1) = sqrt(pg(i)) * w;
        end
        t = (0:Nfft-1)*Ts;
    end
    
    function [tau_mean, tau_rms] = rmsDelaySpread(pg, tau)

        % Normalize PDP (important)
        pg = pg / sum(pg);
    
        % Mean delay
        tau_mean = sum(pg .* tau);
    
        % RMS delay spread
        tau_rms = sqrt(sum(pg .* (tau - tau_mean).^2));
    
    end
end

end