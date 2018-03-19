function [acqResults,result] = acquisition1(rawSignal,longSignal, fftCode, settings, interTime)
%% Initialization =========================================================
% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
% Create two 1msec vectors of data to correlate with and one with zero DC
% signal = (longSignal(1:12000*samplesPerCode).*(longSignal(25:12000*samplesPerCode+24)));
rawSignal = double(rawSignal(1:101*samplesPerCode));
satNum = length(settings.acqSatelliteList);
%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq   = zeros(1, satNum);
% C/A code phases of detected signals
acqResults.codePhase  = zeros(1, satNum);
% Correlation peak ratios of the detected signals
acqResults.peakMetric = zeros(1, satNum);
acqResults.codeDopple = zeros(1, satNum);
% intergral time
% interTime = 12000;
msOfTimeIndex = (round((51.15)./[-6.5:0.5:6.5]));
msOfTimeIndex((abs(msOfTimeIndex)>interTime)&(~isinf(msOfTimeIndex))) = [];
dopplerCode = (51.15)./msOfTimeIndex;
result  =  zeros(length(settings.acqSatelliteList),length(msOfTimeIndex),samplesPerCode);
fprintf('(');
    %% Correlate signals ======================================================
for i = 1:length(msOfTimeIndex)
    if isinf(msOfTimeIndex(i))
        corrSignal = sum(reshape(longSignal(1:interTime*samplesPerCode),samplesPerCode,interTime),2)';
    elseif msOfTimeIndex(i) < 0
        m = - msOfTimeIndex(i);
        corrTime   = abs(round(interTime/(2 * m)));
        lenOfTime  = 2 * m * samplesPerCode + 1;
        corrSignal = sum(reshape(longSignal(1:corrTime*lenOfTime),lenOfTime,corrTime),2);
        corrSignal(m * samplesPerCode + 1) = [];
        corrSignal = sum(reshape(corrSignal,samplesPerCode,2*m),2)';
    else
        m = msOfTimeIndex(i);
        corrTime   = round(interTime/(2*m));
        lenOfTime  = 2 * m * samplesPerCode - 1;
        corrSignal = sum(reshape(longSignal(1:corrTime * lenOfTime), lenOfTime, corrTime),2);
        corrSignal = sum(reshape(corrSignal([1:m*samplesPerCode m*samplesPerCode:end]),samplesPerCode,2*m),2)';
    end
    signalFreqDom = fft(corrSignal);
    result(:,i,:) = bsxfun(@times,signalFreqDom,fftCode);
end
result = abs(ifft(result,samplesPerCode,3)).^2;    
%% Look for correlation peaks in the results ==============================
[~ ,codeDoppleIndex ] = max(max(result,[],3),[],2);
[peakSize ,codePhase] = max(max(result,[],2),[],3); 
acqResults.codePhase  = codePhase;
acqResults.codeDopple = dopplerCode(codeDoppleIndex);
%--- Find 1 chip wide C/A code phase exclude range around the peak ----
samplesPerCodeChip = round(settings.samplingFreq / settings.codeFreqBasis);
for ii = 1:satNum 
    id1 = codePhase(ii) - samplesPerCodeChip;
    id2 = codePhase(ii) + samplesPerCodeChip;
    A   = 1:samplesPerCode;
    codePhaseRange = A((A<id1)&(A>id2-samplesPerCode)|((A>id2)&(A<id1+samplesPerCode)));
    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(result(ii,codeDoppleIndex(ii),codePhaseRange));    
    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(ii) = peakSize(ii)/secondPeakSize;
    % If the result is above threshold, then there is a signal ...
    if  acqResults.peakMetric(ii) > 2.2
        %% Fine resolution frequency search =======================================
        %--- Indicate PRN number of the detected signal -------------------
        PRN = settings.acqSatelliteList(ii);
        fprintf('%02d ', PRN);
        caCodes  = makeCaTable(settings,30,acqResults.codeDopple(ii),PRN);
        xCarrier = rawSignal(codePhase(ii):codePhase(ii) + 30*samplesPerCode - 1) .* caCodes;
        xCarrier = reshape(xCarrier,samplesPerCode,30);
        cztxc    = abs(czt(xCarrier,201,exp(-1i*2*pi*10/settings.samplingFreq),...
            exp(1i*2*pi*(acqResults.codeDopple(ii)*763-1000+settings.IF)/settings.samplingFreq)));
        cztxc    = sum(cztxc,2);
        cztFreqBins = (acqResults.codeDopple(ii)*763-1000+settings.IF):10:...
            (acqResults.codeDopple(ii)*763+1000+settings.IF);
        [~, cztMaxIndex] = max(cztxc);
        figure();
        plot(cztFreqBins,cztxc);
        %--- Save properties of the detected satellite signal -------------
        acqResults.carrFreq(ii)  = cztFreqBins(cztMaxIndex);      
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end
end% for PRN = satelliteList
%=== Acquisition is over ==================================================
fprintf(')\n');
