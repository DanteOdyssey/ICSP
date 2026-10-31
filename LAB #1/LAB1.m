clear; close all; clc;

%% SPECS
fp = 1e6;             % Pass-band edge [Hz]
BW = 3e6;             % Transition-band width [Hz]
fs = fp + BW;         % Stop-band edge [Hz]
N = 4;                % Filter Order
Rp = 0.5;             % Maximum pass-band ripple [dB]
Rs = 50;              % Minimum stop-band attenuation [dB]

%% ELLIPTIC FILTER DESIGN

[b,a] = ellip(N, Rp, Rs, 2*pi*fp, 's');
H = tf(b,a);

%% MAGNITUDE PLOT

f = logspace(3, 8, 10000);  % Frequency [Hz]
w = 2*pi*f;                 % Angular frequency [rad/s]

Hjw = freqs(b, a, w);
Mag_dB = 20*log10(abs(Hjw));

figure;
semilogx(f, Mag_dB, 'b', 'LineWidth', 1.5);
grid on;

xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
title('4th-Order Elliptic Low-Pass Filter');

xlim([min(f) max(f)]);
ylim([-80 5]);

figure;
semilogx(f, Mag_dB, 'b', 'LineWidth', 1.5);
grid on;

xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
title('4th-Order Elliptic Low-Pass Filter - Passband Zoom');

xlim([min(f) fp]);
ylim([-0.6 0.1]);

figure;
semilogx(f, Mag_dB, 'b', 'LineWidth', 1.5);
grid on;

xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
title('4th-Order Elliptic Low-Pass Filter - Rejection Zone Zoom');

xlim([fs max(f)]);
ylim([-100 -50]);

%% FIND BIQUAD COEFFICIENTS
[z,p,k] = tf2zp(b,a);
sos = zp2sos(z, p, k);

fprintf('\n===== BIQUAD COEFFICIENTS =====\n');

k0 = zeros(1,2);
k1 = zeros(1,2);
k2 = zeros(1,2);
w0 = zeros(1,2);
Q = zeros(1,2);
for i = 1:2
    fprintf('\nBIQUAD %d\n',i);
    fprintf('Numerator:   %.8g  %.8g  %.8g\n', ...
        sos(i,1),sos(i,2),sos(i,3));
    fprintf('Denominator: %.8g  %.8g  %.8g\n', ...
        sos(i,4),sos(i,5),sos(i,6));
    k0(i) = sos(i,1);
    k1(i) = sos(i,2);
    k2(i) = sos(i,3);
    w0(i) = sqrt(sos(i,6));
    Q(i) = 1 / (w0(i) * sos(i,5));
    fprintf(' k0: %.8g \n k1: %.8g \n k2: %.8g \n w0: %.8g \n Q: %.8g\n', ...
         k0, k1, k2, w0, Q);
end

