%time domain pulse for raised cosine, together with time vector to plot it against
%oversampling factor= how much faster than the symbol rate we sample at
%len=where to truncate response (multiple of symbol time) on each side of peak
%a = excess bandwidth
function [rc,time_axis] = raised_cosine(a,m,len)
len_os = floor(len*m); %number of samples on each side of peak
%time vector (in units of symbol interval) on one side of the peak
z = cumsum(ones(len_os,1))/m;
A= sin(pi*(1-a)*z)./z; %term 1
B= cos(pi*(1+a)*z)*4*a; %term 2
C= (1 - (4*a*z).^2)*pi; %term 
zerotest = m/(4*a); %location of zero in denominator
%check whether any sample coincides with zero location
if (zerotest == floor(zerotest))
% for taking care of undefined values
B(zerotest) = cos(pi*(1+a)*(z(zerotest)+0.001))*4*a;
A(zerotest) = sin(pi*(1-a)*(z(zerotest)+0.001))./(z(zerotest)+0.001); %term 2
C(zerotest) = (1-(4*a*(z(zerotest)+0.001)).^2)*pi;

end
D = (A+B)./C; %response to one side of peak
rc = [flipud(D);1;D]; %add in peak and other side
rc(len_os+1)=(4*a/pi) + (1-a);
time_axis = [flipud(-z);0;z];
end