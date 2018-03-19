function output_int16
file = 'E:\data\20171114_4M_indoor.bin';
outfile  = 'E:\data\20171114_raw_4M_indoor(2).bin';
fidread  = fopen(file,'r');
fidwrite = fopen(outfile,'w');
N = 110;
numPerfile = 20000000;
I  = kron(ones(1,numPerfile/8,'int8'),int8([1,0,-1,0]));
Q  = kron(ones(1,numPerfile/8,'int8'),int8([0,1,0,-1]));
hd = bdsfilter_4M_2;
for i=1:N
    data  = fread(fidread, numPerfile,'int16=>int8');
    data  = reshape(data, 2, numPerfile/2);
    rawY  = data(1,:).*I - data(2,:).*Q;
    y     = filter(hd,[rawY zeros(1,34)]);
    output_mark = fwrite(fidwrite, y(35:end), 'int8');
    if(output_mark<0)
        disp('output failure!');
        break;
    end
end
fclose('all');    