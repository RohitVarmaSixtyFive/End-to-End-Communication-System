
%%
clear vars;
clc;

[~, fs] = audioread('project.wav');
wavdata = [120/128,250/128,10/128];

snr = 10;

out = zeros(1,length(snr));
for z = 1:length(snr)
%% Quantization
range = 2; 
level8 = range / (2^8);
eight = round(wavdata / level8) * level8;

figure;
subplot(4,1,1);
stem(eight);
Tb = 1/(fs*9);
%% BPSK conversion
encoded = bpsk_map(eight);
subplot(4,1,2);
stem(encoded);

%% line coding 
%rasied cosine
a = 1; 
m = 9; 
len = 2;
[rc,time] = raised_cosine(a,m,len);
rc = rc.*(1/max(rc));

encoded_upsample = upsample(encoded,length(rc));
% rectanngular pulse
for i=1:length(rc)
    if(i>=1&&i<=37)
    rc(i) = 1;
    else
    rc(i) = 0;
    end
end
%%
out_line = conv(rc,encoded_upsample);
out_line = out_line(1:(length(out_line)-(length(rc)-2)));
subplot(x,1,3);
plot(out_line);

%% Modulation
fc = 1e6;
t = 0:1/(10*fc):(length(out_line)-1)/(10*fc);
modulated = out_line.*cos(2*pi*fc*t)';
subplot(x,1,4);
plot(cos(2*pi*fc*t));
subplot(x,1,5);
plot(modulated);

%% noise
EbNo = snr(z);
EbNo_linear = 10^(EbNo / 10);

gt_temp = cumsum(modulated(1:length(rc)).*modulated(1:length(rc)));
gt_norm = gt_temp(length(gt_temp));
len_bit = 1;
Eb = gt_norm*6/len_bit/length(rc);

No = Eb/EbNo_linear;
sigma = sqrt(No/2);

% noisy_modulated = modulated + sigma*randn(length(modulated),1);
% subplot(x,1,6); 
% plot(noisy_modulated);

%% memory noise

delta = zeros(length(rc)+1,1);
a = 0.8;
delta(1) = a;
delta(end) = 1-a;
modulated = conv(modulated,delta);
modulated = modulated(1:end-length(rc));

noisy_modulated = modulated + sigma*randn(length(modulated),1);
subplot(x,1,6); 
plot(noisy_modulated);

%% demodulation

demod = 2*cos(2*pi*fc*t);
demodulated_blp = noisy_modulated.*demod';
% demodulated_blp = modulated.*demod';
demodulated_lp = lowpass(demodulated_blp,1000,fs);
% demodulated_lp = demodulated_blp;

figure;
subplot(3,1,1);
plot(demodulated_lp);

%% line decode
demodulated_lp_cluster = reshape(demodulated_lp(1:end-1),length(rc),[]);

line_decoded = zeros(size(demodulated_lp_cluster, 2), 1);

for i = 1:size(demodulated_lp_cluster, 2)
    row = demodulated_lp_cluster(:, i);
    [maximum, idx] = max(abs(row));

    line_decoded(i) = row(19);
end
subplot(3,1,2);
stem(line_decoded);

%% decoder
decoded = zeros(1,length(line_decoded));
for i = 1:length(decoded)
    if(line_decoded(i)>0)
        decoded(i) = 1;
    else
        decoded(i) = 0;
    end
end

subplot(3,1,3);
stem(decoded);

%% digital to analog
decoded_temp = reshape(decoded, 9, []);
da_out = zeros(1, size(decoded_temp, 2));

for i = 1:size(decoded_temp, 2)
    binary_string = num2str(decoded_temp(:, i)');
    binary_string = strrep(binary_string, ' ', '');  % to remove spaces ig,might not be required for every input but lite
    da_out(i) = bin2dec(binary_string)/128;
end
% Calculate deviations
deviation = 0;
for i = 1:length(da_out)
    if(abs(da_out(i) - eight(i)) > 1e-6)
        deviation = deviation + 1;
        index = i;
    end
end
disp(['Total deviations: ', num2str(deviation/length(da_out))]);
da_out = da_out-1;
out(z) = deviation/length(da_out);

end

