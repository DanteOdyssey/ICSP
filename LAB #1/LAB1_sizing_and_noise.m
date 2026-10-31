clear; close all; clc;

%% 1) SPECS
fp = 1e6;             % Pass-band edge [Hz]
BW = 3e6;             % Transition-band width [Hz]
fs = fp + BW;         % Stop-band edge [Hz]
N  = 4;               % Filter order
Rp = 0.5;             % Maximum pass-band ripple [dB]
Rs = 50;              % Minimum stop-band attenuation [dB]

%% 2) CHOICES
% capacitors, one per biquad in cascade order
Ca = [15e-12 15e-12];     % OP2 [F]
Cb = [15e-12 15e-12];     % OP1 [F]

% true: low-Q biquad first
lowQfirst = true;

% max value of one RNHR1000 resistor [Ohm], Inf = no check
Rmax = Inf;

%% 3) ELLIPTIC FILTER AND POLE/ZERO PAIRS
[b,a] = ellip(N, Rp, Rs, 2*pi*fp, 's');
[z,p,k] = tf2zp(b,a);
Hdc = abs(polyval(b,0)/polyval(a,0));      % DC gain of the prototype (-Rp dB)

pp = p(imag(p) > 0);   pp = pp(:).';       % one pole per complex-conjugate pair
zz = z(imag(z) > 0);   zz = zz(:).';       % one zero per pair (jwz)
nst = numel(pp);

w0_all = abs(pp);
Q_all  = w0_all ./ (-2*real(pp));          % Q = w0/(2|Re p|)
[Q_all, ix] = sort(Q_all, 'ascend');       % index 1 = lowest Q
w0_all = w0_all(ix);
wz_all = sort(abs(zz), 'ascend');          % index 1 = zero closest to the pass band
wz_pair = fliplr(wz_all);                  % high Q <-> zero closest to the pass band

if lowQfirst
    ord = 1:nst;
else
    ord = nst:-1:1;
end
w0 = w0_all(ord);
Q  = Q_all(ord);
wz = wz_pair(ord);

g = ones(1,nst);  g(end) = Hdc;    % prototype DC gain on the last stage
k0 = g;                            % DC gain
k1 = zeros(1,nst);                 % no s term in the numerator
k2 = g ./ wz.^2;                   % k2 = k0/wz^2  [s^2]

%% 4) SIZING: FROM C TO R
R1 = Q ./ (w0 .* Ca);
R3 = 1 ./ (w0 .* Ca);
R4 = 1 ./ (w0 .* Cb .* k0);
R5 = -1 ./ (w0 .* Cb);
C4 = -(w0.^2) .* k2 .* Ca;
R2 = inf(1,nst);                   % k1 = 0 -> branch not present

fprintf('BIQUADS (in cascade order):\n');
for i = 1:nst
    fprintf('Biquad %d: f0 = %.4f MHz | Q = %.4f | fz = %.4f MHz | k0 = %.5f\n', i, w0(i)/2/pi/1e6, Q(i), wz(i)/2/pi/1e6, k0(i));
end

crossStr = {'', '   [swap wires]'};
fprintf('\nCOMPONENTS (absolute values):\n');
for i = 1:nst
    fprintf('\nBIQUAD %d   (Ca = %.4g pF, Cb = %.4g pF)\n', i, Ca(i)*1e12, Cb(i)*1e12);
    fprintf('  R1 = %10.4f kOhm%s\n', abs(R1(i))/1e3, crossStr{(R1(i) < 0) + 1});
    fprintf('  R3 = %10.4f kOhm%s\n', abs(R3(i))/1e3, crossStr{(R3(i) < 0) + 1});
    fprintf('  R4 = %10.4f kOhm%s\n', abs(R4(i))/1e3, crossStr{(R4(i) < 0) + 1});
    fprintf('  R5 = %10.4f kOhm%s\n', abs(R5(i))/1e3, crossStr{(R5(i) < 0) + 1});
    fprintf('  C4 = %10.4f pF  %s\n',  abs(C4(i))*1e12, crossStr{(C4(i) < 0) + 1});
end

if isfinite(Rmax)                  % resistors to be put in series
    Rall = abs([R1; R3; R4; R5]);
    names = {'R1','R3','R4','R5'};
    for i = 1:nst
        for r = 1:4
            n = ceil(Rall(r,i)/Rmax);
            if n > 1
                fprintf('Biquad %d, %s = %.4f kOhm -> %d x %.4f kOhm in series\n',i, names{r}, Rall(r,i)/1e3, n, Rall(r,i)/n/1e3);
            end
        end
    end
end

%% 5) CHECK: H REBUILT FROM THE COMPONENTS vs IDEAL ELLIPTIC
f = logspace(3, 8, 5000);
s = 1j*2*pi*f;
Hst = zeros(nst, numel(f));
for i = 1:nst
    num = -(C4(i)/Ca(i))*s.^2 - s*(1/R2(i))/Ca(i) + 1/(R3(i)*R4(i)*Ca(i)*Cb(i));
    den = s.^2 + s/(R1(i)*Ca(i)) - 1/(R3(i)*R5(i)*Ca(i)*Cb(i));
    Hst(i,:) = num ./ den;
end
Hchain = prod(Hst, 1);
Href = freqs(b, a, 2*pi*f);
err_dB = max(abs(20*log10(abs(Hchain)) - 20*log10(abs(Href))));
fprintf('\nMax error |H_components| vs |H_ellip| = %.2e dB (should be ~0)\n', err_dB);

bq = 'gm';
figure;
semilogx(f, 20*log10(abs(Href)), 'c', 'LineWidth', 1.5); hold on;
semilogx(f, 20*log10(abs(Hchain)), 'y--', 'LineWidth', 1.2);
for i = 1:nst
    semilogx(f, 20*log10(abs(Hst(i,:))), 'Color', bq(mod(i-1,numel(bq))+1), 'LineWidth', 1);
end
grid on; xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]'); ylim([-100 20]);
legend([{'Ideal ellip', 'From components'}, ...
        arrayfun(@(i) sprintf('Biquad %d', i), 1:nst, 'UniformOutput', false)],'Location', 'southwest');
title('Biquad and filter responses');

%% 6) MASK CHECK (red = spec limits)
Mag_dB = 20*log10(abs(Href));       % ideal ellip
Mag_c  = 20*log10(abs(Hchain));     % from components
lim = @(y) plot([f(1) f(end)], [y y], 'r--', 'LineWidth', 1);

figure;                              % full response
semilogx(f, Mag_dB, 'c', 'LineWidth', 1.5); hold on;
semilogx(f, Mag_c, 'y--', 'LineWidth', 1);
lim(-Rp); lim(-Rs);
xline(fp, 'r:'); xline(fs, 'r:');
grid on; xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
xlim([min(f) max(f)]); ylim([-100 5]);
legend('Ideal ellip', 'From components', 'Spec limits', 'Location', 'southwest');
title('Complete filter and mask');

figure;                              % pass-band zoom
semilogx(f, Mag_dB, 'c', 'LineWidth', 1.5); hold on;
semilogx(f, Mag_c, 'y--', 'LineWidth', 1);
lim(0); lim(-Rp);
grid on; xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
xlim([min(f) fp]); ylim([-0.6 0.1]);
title('Pass-band zoom (ripple between 0 and -Rp)');

figure;                              % stop-band zoom
semilogx(f, Mag_dB, 'c', 'LineWidth', 1.5); hold on;
semilogx(f, Mag_c, 'y--', 'LineWidth', 1);
lim(-Rs);
grid on; xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
xlim([fs max(f)]); ylim([-100 -40]);
title('Stop-band zoom (below -Rs)');





%% 7) NOISE TRANSFER FUNCTIONS (one set per biquad)
fn = logspace(4, 9, 10000);            % own frequency axis (slides: 1e5 ... 1e9)
sn = 1j*2*pi*fn;

% Physical resistor values
r1 = abs(R1);  r3 = abs(R3);  r4 = abs(R4);  r5 = abs(R5);

HnR1 = zeros(nst, numel(fn));   HnR3 = HnR1;   HnR4 = HnR1;
HnR5 = HnR1;                    HnA1 = HnR1;   HnA2 = HnR1;

for i = 1:nst
    % Common denominator: s^2/w0^2 + s/(w0*Q) + 1
    D = sn.^2/w0(i)^2 + sn/(w0(i)*Q(i)) + 1;

    par45 = r4(i)*r5(i) / (r4(i) + r5(i));     % R4 // R5
    par13 = r1(i)*r3(i) / (r1(i) + r3(i));     % R1 // R3

    HnR1(i,:) = (sn/(w0(i)*Q(i)))              ./ D;
    HnR3(i,:) = (sn/w0(i) * Cb(i)/Ca(i))       ./ D;
    HnR4(i,:) = k0(i)                          ./ D;
    HnR5(i,:) = 1                              ./ D;
    HnA1(i,:) = ((r4(i)+r5(i))/r4(i) * (1 + sn*Ca(i)*par45)) ./ D;
    HnA2(i,:) = ((r1(i)+r3(i))/r1(i) * r5(i) * Ca(i) * sn ...
        .* (1 + sn*Cb(i)*par13)) ./ D;
end

%% 8) PLOTS 
names = {'R1 noise', 'R3 noise', 'R4 noise', 'R5 noise', ...
    'First amplifier noise', 'Second amplifier noise'};
for i = 1:nst
    Hn = {HnR1(i,:), HnR3(i,:), HnR4(i,:), HnR5(i,:), HnA1(i,:), HnA2(i,:)};

    figure('Name', sprintf('Noise TFs - biquad %d', i));
    tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Biquad %d  (f0 = %.3f MHz, Q = %.2f)', i, w0(i)/2/pi/1e6, Q(i)));

    for m = 1:6
        nexttile;
        semilogx(fn, 20*log10(abs(Hn{m}) + eps), 'r', 'LineWidth', 1.3);
        grid on; xline(w0(i)/2/pi, 'k:');
        xlabel('Frequency [Hz]'); ylabel('|V_o/V_n| [dB]');
        title(names{m});
    end
end



%% 9) RESISTOR NOISE POWER (integral from 1 Hz to 10 GHz)
kB = 1.380649e-23;   T = 300;                 % [J/K], [K]
toFilterOutput = true;   % false: output of the biquad that contains the resistor
% true : refer to the filter output (multiplies by the following biquads)

fi = logspace(0, 10, 2e5);                    % 1 Hz ... 10 GHz [Hz]
si = 1j*2*pi*fi;                              % s = j*2*pi*f

srcName = {'R1','R3','R4','R5'};
Pres = zeros(nst, 4);                         % integrated power [V^2], rows = biquads

for i = 1:nst
    D = si.^2/w0(i)^2 + si/(w0(i)*Q(i)) + 1;  % common denominator

    % Gain of the following biquads (only if referring to the filter output)
    Gnext = ones(size(fi));
    if toFilterOutput
        for j = i+1:nst
            Gnext = Gnext .* (k0(j) + k2(j)*si.^2) ./ (si.^2/w0(j)^2 + si/(w0(j)*Q(j)) + 1);
        end
    end

    % Noise transfer functions (complex), same formulas as section 7
    Hr = [ (si/(w0(i)*Q(i)))        ./ D ;
        (si/w0(i)*Cb(i)/Ca(i))   ./ D ;
        k0(i)*ones(size(fi))     ./ D ;
        ones(size(fi))           ./ D ] .* Gnext;

    % Voltage noise PSD of the resistors [V^2/Hz]: 4kTR
    Sv = 4*kB*T*[r1(i); r3(i); r4(i); r5(i)];

    % Output PSD per source: |H(j2*pi*f)|^2 * 4kTR  [V^2/Hz]
    Sout = abs(Hr).^2 .* Sv .* 2;                  % 4 x numel(fi)

    % Integration in f (Hz): total and cumulative
    Pres(i,:) = trapz(fi, Sout, 2).';
    Pcum      = cumtrapz(fi, Sout, 2);        % cumulative integral from 1 Hz

    % Cumulative integral plot
    figure('Name', sprintf('Resistor noise integration - biquad %d', i));
    semilogx(fi, Pcum.'*1e12, 'LineWidth', 1.2); hold on;
    semilogx(fi, sum(Pcum,1)*1e12, 'k--', 'LineWidth', 1.2);
    grid on; xlim([1e5 1e9]);
    xlabel('Frequency [Hz]'); ylabel('Cumulative noise power [\muV^2]');
    legend([srcName, {'Total'}], 'Location', 'northwest');
    title(sprintf('Biquad %d  (f0 = %.3f MHz, Q = %.2f): cumulative integral from 1 Hz', ...
        i, w0(i)/2/pi/1e6, Q(i)));
end

%% Results
Ptot_bq  = sum(Pres, 2);          % total resistor noise power of each biquad [V^2]
Ptot_all = sum(Pres(:));          % total over all biquads [V^2]

fprintf('\nRESISTOR NOISE POWER (1 Hz - 10 GHz), T = %g K', T);
if toFilterOutput, fprintf('  [referred to filter output]\n');
else,              fprintf('  [at each biquad output]\n'); end

for i = 1:nst
    fprintf('\nBiquad %d  (total = %.4e V^2 -> %.2f uV rms, %.2f %% of overall total)\n', ...
        i, Ptot_bq(i), sqrt(Ptot_bq(i))*1e6, 100*Ptot_bq(i)/Ptot_all);
    for m = 1:4
        fprintf('  %s: %.4e V^2 (%7.2f uV rms) | %6.2f %% of biquad | %6.2f %% of overall total\n', ...
            srcName{m}, Pres(i,m), sqrt(Pres(i,m))*1e6, ...
            100*Pres(i,m)/Ptot_bq(i), 100*Pres(i,m)/Ptot_all);
    end
end

fprintf('\nOverall total (all biquads, all resistors) = %.4e V^2 -> %.2f uV rms\n', ...
    Ptot_all, sqrt(Ptot_all)*1e6);