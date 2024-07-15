%% Clearing variables and console
clear vars;
clc;
%% Test Input for all the plots
[~, fs] = audioread('project.wav');
wavdata = [120/128, 250/128, 10/128, 115/128];

%% SNR Configuration
snr = 1000;

%% SNR Output
out = zeros(1, length(snr));

%% Main SNR loop
for z = 1:length(snr)
    %% Quantization 
    range = 2;
    level8 = range / (2^8);
    eight = round(wavdata / level8) * level8;
    Tb = 1/(fs*9);

    %% A/D Conversion
    scaled_in = round(eight * 128);
    ad_out = zeros(length(scaled_in) *9,1);

    for i = 1:length(scaled_in)
        binStr = dec2bin(scaled_in(i), 9);
        binArray = binStr - '0'; 
        ad_out((i-1)*9 + 1:i*9) = binArray; 
    end
    encoded = ad_out;

    %% BPSK Conversion
    for i = 1:length(encoded)
        if(encoded(i)==0)
            encoded(i) = -1;
        end
    end

    %% Line Coding - Raised Cosine
    a = 1;
    m = 9;
    len = 2;
    [rc, time] = raised_cosine(a, m, len);
    rc = rc.*(1/max(rc));
    encoded_upsample = upsample(encoded, length(rc));

    %% Rectangular Pulse
    for i = 1:length(rc)
        if (i >= 1 && i <= 37)
            rc(i) = 1;
        else
            rc(i) = 0;
        end
    end

    %% Line Output
    out_line = conv(rc, encoded_upsample);
    out_line = out_line(1:(length(out_line) - (length(rc) - 2)));

    %% Modulation
    fc = 1e6;
    t = 0:1/(10*fc):(length(out_line) - 1)/(10*fc);
    modulated = out_line.*cos(2*pi*fc*t)';

    %% Noise Addition
    EbNo = snr(z);
    EbNo_linear = 10^(EbNo / 10);
    gt_temp = cumsum(modulated(1:length(rc)).*modulated(1:length(rc)));
    gt_norm = gt_temp(length(gt_temp));
    len_bit = 1;
    Eb = gt_norm*6/len_bit/length(rc);
    No = Eb/EbNo_linear;
    sigma = sqrt(No/2);
    noisy_modulated = modulated + sigma*randn(length(modulated), 1);

    %% Memory Noise
    delta = zeros(length(rc) + 1, 1);
    a = 0.8;
    delta(1) = a;
    delta(end) = 1-a;
    modulated = conv(modulated, delta);
    modulated = modulated(1:end-length(rc));

    %% Demodulation
    cutoff = 1000;
    demod = 2*cos(2*pi*fc*t);
    demodulated_blp = noisy_modulated.*demod';
    demodulated_lp = lowpass(demodulated_blp, cutoff, fs);

    %% Line Decode
    demodulated_lp_cluster = reshape(demodulated_lp(1:end-1), length(rc), []);
    line_decoded = zeros(size(demodulated_lp_cluster, 2), 1);
    for i = 1:size(demodulated_lp_cluster, 2)
        row = demodulated_lp_cluster(:, i);
        line_decoded(i) = row(19);
    end

    %% Decoder
    decoded = zeros(1, length(line_decoded));
    for i = 1:length(decoded)
        if (line_decoded(i) > 0)
            decoded(i) = 1;
        else
            decoded(i) = 0;
        end
    end

    %% Digital to Analog Conversion
    decoded_temp = reshape(decoded, 9, []);
    da_out = zeros(1, size(decoded_temp, 2));
    for i = 1:size(decoded_temp, 2)
        binary_string = num2str(decoded_temp(:, i)');
        binary_string = strrep(binary_string, ' ', '');
        da_out(i) = bin2dec(binary_string) / 128;
    end

    %% Calculate Deviations
    deviation = 0;
    for i = 1:length(da_out)
        if (abs(da_out(i) - eight(i)) > 1e-6)
            deviation = deviation + 1;
        end
    end
    disp(['Total deviations: ', num2str(deviation / length(da_out))]);
    da_out = da_out - 1;
    out(z) = deviation / length(da_out);
end

%% All Plots
    figure;
    subplot(4,1,1);
    stem(ad_out);
    title('A/D Output');
    xlabel('Sample Index');
    ylabel('x_1');
    grid on;
    
    subplot(4,1,2);
    stem(encoded);
    title('BPSK Encoded Signal');
    xlabel('Sample Index');
    ylabel('x_2');
    grid on;
    
    subplot(4,1,3);
    plot(out_line);
    title('Line Coded Signal with Raised Cosine');
    xlabel('Time (s)');
    ylabel('x_3');
    grid on;
    
    subplot(4,1,4);
    plot(modulated);
    title('Modulated Signal');
    xlabel('Time (s)');
    ylabel('x_4');
    grid on;
    
    figure;
    subplot(4,1,1);
    plot(noisy_modulated);
    title('Noisy Modulated Signal');
    xlabel('Time (s)');
    ylabel('y_4');
    grid on;
    
    subplot(4,1,2);
    plot(demodulated_lp);
    title('Demodulated Signal after Lowpass');
    xlabel('Time (s)');
    ylabel('y_3');
    grid on;
    
    subplot(4,1,3);
    stem(line_decoded);
    title('Line Decoded Signal');
    xlabel('Sample Index');
    ylabel('y_2');
    grid on;
    
    subplot(4,1,4);
    stem(decoded);
    title('Binary Decoded Signal');
    xlabel('Sample Index');
    ylabel('y_1');
    grid on;

