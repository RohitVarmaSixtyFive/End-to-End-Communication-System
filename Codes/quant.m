
%%
fs = 400;
f = 25;
Ts = 1/fs;
t = 1:Ts:2;

x = 8*cos(2*pi*f*t);

range = 7; 

level2 = (2*range)/(2^2);
level6 = (2*range)/(2^6);
two = round(x/level2)*level2;
six = round(x/level6)*level6;

qe2 = x - two;
qe6 = x - six;
quantization_error2 = mean(abs(qe2));
disp(['2 bit quantization error: ', num2str(quantization_error2)]);

quantization_error6 = mean(abs(qe6));
disp(['6 bit quantization error: ', num2str(quantization_error6)]);
%% b 
figure();
subplot(3,2,1);
plot(t,x);

subplot(3,2,3);
plot(t,six);
title('Plot for 6 bit quantization')

subplot(3,2,5);
plot(t,qe6);
title('Plot for 6 bit quantization error')

%% b
subplot(3,2,2);
plot(t,x);
title('Plot for x')

subplot(3,2,4);
plot(t,two);
title('Plot for 2 bit quantization')

subplot(3,2,6);
plot(t,qe2);
title('Plot for 2 bit quantization error')

%% c
snr2 = qe2.*qe2;
snr6 = qe6.*qe6;

Signal_energy = x.*x;

final_snr2 = mean(Signal_energy)/mean(snr2);
final_snr6 = mean(Signal_energy)/mean(snr6);

db_final_snr2 = 10*log10(final_snr2);
db_final_snr6 = 10*log10(final_snr6);

disp(['2 bit snr(db): ', num2str(db_final_snr2)]);
disp(['6 bit snr(db): ', num2str(db_final_snr6)]);
disp(['2 bit snr: ', num2str(final_snr2)]);
disp(['6 bit snr: ', num2str(final_snr6)]);

