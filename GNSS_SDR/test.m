function  s = test( x )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
a=[-3.2286,-1.4796,-2.0482,-3.9590;2.7106,3.9479,2.1879,1.4503;0.11345,0.070453,2.9621,0.05787]*10^7;
s(1)=pdist2(x(1:3),a(:,1)')-x(4);
s(2)=pdist2(x(1:3),a(:,2)')-x(4);
s(3)=pdist2(x(1:3),a(:,3)')-x(4);
s(4)=pdist2(x(1:3),a(:,4)')-x(4);

end
