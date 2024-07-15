clear vars;
clc;
[wavdata, fs] = audioread('project.wav');
wavdata = wavdata(:,1);
wavdata = wavdata + 1; % shift lil up to make values positive. as bin2dec is doubling bit length for neg values
% [~, fs] = audioread('project.wav');
% wavdata = [120/128,250/128,10/128,115/128];
wavdata = wavdata(1:100000);

%% Quantization
range = 2; 
level8 = range / (2^8);
eight = round(wavdata / level8) * level8;
% eight = [150/128; 180/128];
x = 6;
Tb = 1/(9*fs);
%% BPSK conversion
encoded = bpsk_map(eight);
%% line coding 

%rasied cosine
a = 1; 
m = 9; 
len = 2;
[rc,time] = raised_cosine(a,m,len);
rc = rc.*(1/max(rc));

encoded_upsample = upsample(encoded,length(rc));
% rectanngular pulse
% for i=1:length(rc)
%     if(i>=1&&i<=37)
%     rc(i) = 1;
%     else
%     rc(i) = 0;
%     end
% end
%%
out_line = conv(rc,encoded_upsample);
out_line = out_line(1:(length(out_line)-(length(rc)-2)));


%% Modulation
fc = 1e6;
t = 0:1/(10*fc):(length(out_line)-1)/(10*fc);
modulated = out_line.*cos(2*pi*fc*t)';

%% noise
EbNo = 0;
EbNo_linear = 10^(EbNo / 10);

gt_temp = cumsum(modulated(1:length(rc)).*modulated(1:length(rc)));
gt_norm = gt_temp(length(gt_temp));
Eb = gt_norm/length(rc);

No = Eb/EbNo_linear;
sigma = sqrt(No/2);

sigma = 0.4;
% noisy_modulated = modulated + sigma*randn(length(modulated),1) + 1j*sigma*randn(length(modulated),1);
 %% Memory Noise
b = 2;
    delta = zeros(b*length(rc) + 1, 1);
    a = 0.8;
    delta(1) = a;
    delta(end) = 1-a;
    modulated = conv(modulated, delta);
    modulated = modulated(1:end-length(rc));
    noisy_modulated_mem = modulated + sigma*randn(length(modulated), 1) + 1j*sigma*randn(length(modulated), 1);

%% demodulation

demod = 2*cos(2*pi*fc*t);
demodulated_blp = noisy_modulated.*demod';
% demodulated_blp = modulated.*demod';
cutoff = fc+1000;  
% order = 6;  
% normalized_cutoff = min(cutoff/(fs/2), 0.999);
% [b, a] = butter(order, normalized_cutoff, 'low');
% demodulated_lp = filter(b, a, demodulated_blp);
demodulated_lp = lowpass(demodulated_blp,1000,fs);



%% line decode
matched_output = conv(rc,demodulated_lp);

t = 0:1/(10*fc):(length(matched_output)-1)/(10*fc);

demodulated_lp_cluster = reshape(matched_output(19:end-19),length(rc),[]);

line_decoded = zeros(size(demodulated_lp_cluster, 2), 1);

for i = 1:size(demodulated_lp_cluster, 2)
    row = demodulated_lp_cluster(:, i);
    [maximum, idx] = max(abs(row));
    
    line_decoded(i) = row(19);
end

    demodulated_lp_cluster_mem = reshape(demodulated_lp_mem(1:end-1), length(rc), []);
    line_decoded_mem = zeros(size(demodulated_lp_cluster_mem, 2), 1);
    for i = 1:size(demodulated_lp_cluster_mem, 2)
        row = demodulated_lp_cluster_mem(:, i);
        line_decoded_mem(i) = row(19);
    end

    %% Constellation plot
    % figure;
    % % Plot for the line decoded signal in AWGN channel
    % subplot(2,1,1);
    % plot(real(line_decoded), imag(line_decoded), '*');
    % title('Output Constellation Diagram for AWGN Channel');
    % xlabel('In-phase Component (I)');
    % ylabel('Quadrature Component (Q)');
    
    % Plot for the line decoded signal in AWGN channel with memory
    figure;
    plot(real(line_decoded_mem), imag(line_decoded_mem), '*');
    title('Output Constellation Diagram for AWGN Channel with Memory with a = 0.8 and b=2');
    xlabel('In-phase Component (I)');
    ylabel('Quadrature Component (Q)');

%% decoder

decoded = zeros(1,length(line_decoded));
for i = 1:length(decoded)
    if(line_decoded(i)>0)
        decoded(i) = 1;
    else
        decoded(i) = 0;
    end
end

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



