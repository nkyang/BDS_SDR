function [ result ] = deinter( data )
%对电文每个字（30bits）解交织
%   [ result ] = deinter( data )
%  
index   = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30];
temp = inf(1,30);
for i = 1:30
    temp(i)=data(index(i));
end
result1 = BCH (temp( 1:15));
result2 = BCH (temp(16:30));
result  = [result1(1:11) result2(1:11) result1(12:15) result2(12:15)];