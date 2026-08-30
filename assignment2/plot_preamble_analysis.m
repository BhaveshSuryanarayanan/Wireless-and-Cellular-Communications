%% Preamble Symbol Analysis
% (a) Frequency-domain and time-domain properties of the preamble

function plot_preamble_analysis()
%%%%% OFDM parameters %%%%%
params.del_f = 1e4;              % subcarrier spacing = 10kHz
params.numCarriers = 512;        % fft size = 512
params.prefixDuration = 6.25e-6; % CP duration = 6.25 μs

params.sampleRate = params.del_f * params.numCarriers;
params.prefixLength = round(params.prefixDuration * params.sampleRate);
params.blockDuration = (params.numCarriers + params.prefixLength) / params.sampleRate;
params.preamblePeriodicity = 5;
params.guardTones = [241:256 -255:-241];
params.DCSubcarrier = [0];
params.pilots = false;
params.pilotPeriodicity = 8;

%%%%% Channel Parameters %%%%%
params.pathGains = [-3 0 -1 -4 -9 -17];
params.delays = [0 7 13 18 21 26];
params.CFO = 0;                  % No CFO for preamble visualization
params.SNR = 100;                % High SNR for clean visualization

%%%%% Modulation parameters %%%%%
params.modulationType = 'QAM';
params.modulationOrder = 16;
params.numFrames = 1;

% Generate preamble
N = params.numCarriers;
Lcp = params.prefixLength;

S = OFDM.generateOFDMframe(params);
s = OFDM.modulate(S, params);

% Extract first symbol with CP
s1 = s(1 : N + Lcp);

%% Plot 1: Frequency Domain - Preamble Symbol
figure('Position', [100, 100, 900, 600]);
stem(abs(S(:, 1)), 'filled', 'LineWidth', 1.5);
hold on;
grid on; grid minor;

% Mark guard tones
guard_indices = [params.guardTones] + N/2;
guard_indices = guard_indices(guard_indices > 0 & guard_indices <= N);
stem(guard_indices, abs(S(guard_indices, 1)), 'filled', 'LineWidth', 2, 'Color', 'red', 'DisplayName', 'Guard tones');

% Mark DC subcarrier
dc_idx = N/2 + 1; % DC is at center
stem(dc_idx, abs(S(dc_idx, 1)), 'filled', 'LineWidth', 2, 'Color', 'magenta', 'DisplayName', 'DC subcarrier');

xlabel('Subcarrier Index', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Magnitude', 'FontSize', 12, 'FontWeight', 'bold');
title('Preamble Symbol: Frequency Domain', 'FontSize', 13, 'FontWeight', 'bold');

% Add annotation for structure
% text(N/4, max(abs(S(:, 1))) * 0.9, 'Odd subcarriers (Schmidl-Cox)', ...
    % 'FontSize', 10, 'HorizontalAlignment', 'center', ...
    % 'BackgroundColor', 'yellow', 'EdgeColor', 'black');

legend('Data', 'Guard tones', 'DC subcarrier', 'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11);

%% Plot 2: Time Domain - First OFDM Symbol with CP
figure('Position', [1050, 100, 900, 600]);
time_axis = (0:length(s1)-1) / params.sampleRate * 1e6; % Convert to microseconds

plot(time_axis, abs(s1), 'LineWidth', 2, 'Color', 'blue');
hold on;
grid on; grid minor;

% Define window regions (in samples)
cp_start = 0;
cp_end = Lcp;
useful_p1_start = Lcp + 1;
useful_p1_end = Lcp + N/2;
useful_p2_start = Lcp + N/2 + 1;
useful_p2_end = Lcp + N;

% Convert to time axis (in microseconds)
cp_time_start = cp_start / params.sampleRate * 1e6;
cp_time_end = cp_end / params.sampleRate * 1e6;
p1_time_start = useful_p1_start / params.sampleRate * 1e6;
p1_time_end = useful_p1_end / params.sampleRate * 1e6;
p2_time_start = useful_p2_start / params.sampleRate * 1e6;
p2_time_end = useful_p2_end / params.sampleRate * 1e6;

% Fill regions with rectangle patches
max_val = max(abs(s1)) * 1.2;
ax = gca;
rectangle(ax, 'Position', [cp_time_start, 0, cp_time_end - cp_time_start, max_val], ...
    'FaceColor', 'green', 'FaceAlpha', 0.15, 'EdgeColor', 'green', 'LineWidth', 2.5);
rectangle(ax, 'Position', [p1_time_start, 0, p1_time_end - p1_time_start, max_val], ...
    'FaceColor', '#FFD700', 'FaceAlpha', 0.15, 'EdgeColor', 'red', 'LineWidth', 2.5);
rectangle(ax, 'Position', [p2_time_start, 0, p2_time_end - p2_time_start, max_val], ...
    'FaceColor', 'cyan', 'FaceAlpha', 0.15, 'EdgeColor', 'blue', 'LineWidth', 2.5);

% Add labels inside rectangles
text(cp_time_end/2, max_val*0.5, '1: CP', ...
    'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
text((p1_time_start + p1_time_end)/2, max_val*0.5, '2: Part 1', ...
    'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
text((p2_time_start + p2_time_end)/2, max_val*0.5, '3: Part 2', ...
    'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

xlabel('Time (μs)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('|Amplitude|', 'FontSize', 12, 'FontWeight', 'bold');
title('First OFDM Symbol: Time Domain (Absolute Value)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Signal', 'Location', 'best', 'FontSize', 10);

set(gca, 'FontSize', 11);
xlim([0 time_axis(end)]);
ylim([0 max_val]);

%% Add text box with key properties
annotation('textbox', [0.15 0.02 0.7 0.08], ...
    'String', {['Preamble Properties: Odd subcarriers only (Schmidl-Cox) | FFT Size: ' num2str(N)], ...
               ['Subcarrier Spacing: ' num2str(params.del_f/1e3) ' kHz | CP Duration: ' num2str(Lcp) ' samples (' num2str(Lcp/params.sampleRate*1e6, '%.2f') ' μs)']}, ...
    'FontSize', 10, 'EdgeColor', 'black', 'BackgroundColor', 'lightyellow', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

sgtitle('(a) Preamble Symbol Analysis: Frequency and Time Domain Properties', ...
    'FontSize', 14, 'FontWeight', 'bold');

%% Calculate and display maximum unambiguous CFO
fprintf('=== Preamble CFO Measurement Analysis ===\n\n');
fprintf('FFT Size (N): %d\n', N);
fprintf('Subcarrier Spacing (Δf): %.2f kHz\n', params.del_f/1e3);
fprintf('Preamble Type: Schmidl-Cox (odd subcarriers only)\n');
fprintf('Non-zero subcarriers: %d (every 2nd tone)\n\n', N/2);

% Maximum unambiguous CFO is limited by subcarrier spacing / 2
max_unambiguous_cfo = params.del_f / 2;
fprintf('Maximum Unambiguous CFO: ±%.2f kHz\n', max_unambiguous_cfo/1e3);
fprintf('(Limited by Nyquist on subcarrier spacing)\n\n');

% Time-domain period properties
symbol_duration = N / params.sampleRate;
cp_duration = Lcp / params.sampleRate;
total_duration = (N + Lcp) / params.sampleRate;

fprintf('Time Domain Properties:\n');
fprintf('  Useful symbol duration: %.3f μs\n', symbol_duration * 1e6);
fprintf('  Cyclic prefix duration: %.3f μs\n', cp_duration * 1e6);
fprintf('  Total OFDM symbol duration: %.3f μs\n\n', total_duration * 1e6);

fprintf('Autocorrelation Properties (Schmidl-Cox):\n');
fprintf('  The N/2-spaced repetition in the preamble creates a strong\n');
fprintf('  autocorrelation peak at lag N/2, enabling timing and CFO estimation.\n');
fprintf('  Correlation metric M(d) peaks when aligned with symbol boundary.\n');
end