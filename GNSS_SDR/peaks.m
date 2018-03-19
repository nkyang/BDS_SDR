function [fitResult, gof, output]= peaks( data, lowRange, peakNum )
highRange = ceil(max(data)/0.05)*0.05;
scrsz = get(groot,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)]);
%% Ê±ÐòÍ¼
subplot(2,2,1);
times = (hours(0)):seconds(0.1):(days(4)+hours(10)+minutes(21)+seconds(39.5));
plot(times,data);
axis([0 100 floor(min(data)) ceil(highRange)]);
xlabel('time ');
ylabel('signal intensity /dB');
%% ·Ö²¼Í¼
subplot(2,2,2);
histogram(data);
xlabel('signal intensity /dB');
ylabel('number of times ')
%% ÄâºÏ
[N,edges] = histcounts(data,lowRange:0.05:highRange);
probN     = N./sum(N)./0.05;
[xData, yData] = prepareCurveData(edges(2:end), probN );
switch peakNum
    case 1
        ft = fittype('gauss1');
    case 2
        ft = fittype('gauss2');
    case 3
        ft = fittype('gauss3');
    case 4
        ft = fittype('gauss4');
    otherwise
        disp('error')
end
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.TolFun= 1.e-8;
opts.TolX  = 1.e-8;
[fitResult, gof,output] = fit( xData, yData, ft, opts );
coeff = coeffvalues(fitResult);
coeff = reshape(coeff,3,peakNum);
%% ÄâºÏÍ¼
subplot(2,2,3);
plot(fitResult,'m-', xData, yData);
hold on
for i = 1:peakNum
    fun = @(x)coeff(1,i)*exp(-((x-coeff(2,i))/coeff(3,i))^2);
    fplot(fun,[floor(lowRange), ceil(highRange)],'c--');
end
axis([lowRange highRange 0 ceil(max(yData)/0.05)*0.05]);
legend('PDF of signal intensity ','fit result','Location', 'NorthEast');
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',1);
xlabel('signal intensity');
ylabel('Probability Distribution Function');
text(lowRange+0.2,ceil(max(yData)*20)/20-0.02,['RMSE = ',num2str(gof.rmse)]);
text(lowRange+0.2,ceil(max(yData)*20)/20-0.06,['R square = ',num2str(gof.rsquare)]);
grid on
%% ÄâºÏ²Ð²îÍ¼
subplot(2,2,4);
plot(fitResult, xData, yData, 'residuals');
xlim([lowRange highRange]);
xlabel('signal intensity -dB');
ylabel('residuals');
%fitResults = peakfit([edges(2:end)' probN'],(lowRange+highRange)/2,...
%   highRange-lowRange,peaknum);
end
