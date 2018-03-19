function [ mlon , mlat ] = wgs2mgs( wlon , wlat )
a = 6378245.0;
ee = 0.00669342162296594323;
x = wlon - 105;
y = wlat - 35;
radlat = wlat/180*pi;
magic = 1-ee*sin(radlat)^2;
yret = -100+2*x+3*y+0.2*y^2+0.1*x*y+0.2*sqrt(abs(x))+40/3*sin(6*pi*x)...
    +40/3*sin(2*pi*x)+40/3*sin(pi*y)+80/3*sin(pi*y/3)+...
    320/3*sin(pi*y/12)+640/3*sin(pi*y/30);
xret = 300+x+2*y+0.1*x^2+0.1*x*y+0.1*sqrt(abs(x))+40/3*sin(6*pi*x)+...
    40/3*sin(2*pi*x)+40/3*sin(pi*x)+80/3*sin(pi*x/3)+...
    100*sin(pi*x/12)+200*sin(pi*x/30);
dlat = yret*180/((a * (1 - ee)) / (magic * sqrt(magic)) * pi);
dlon = xret*180/(a / sqrt(magic) * cos(radlat) * pi);
mlon=wlon+dlon;
mlat=wlat+dlat;
end

