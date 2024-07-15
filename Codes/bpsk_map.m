function out = bpsk_map(in)
    scaled_in = round(in * 128);
    out = zeros(length(scaled_in) *9,1);

    for i = 1:length(scaled_in)
        binStr = dec2bin(scaled_in(i), 9);
        binArray = binStr - '0'; 
        out((i-1)*9 + 1:i*9) = binArray; 
    end
    for i = 1:length(out)
        if(out(i)==0)
            out(i) = -1;
        end
    end
end
