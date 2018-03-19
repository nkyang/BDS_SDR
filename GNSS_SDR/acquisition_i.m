function [acqResults,result] = acquisition_i(longSignal, fftCode,settings)


%% Initialization =========================================================

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));

% Create two 1msec vectors of data to correlate with and one with zero DC
signal = (longSignal(1:1200*samplesPerCode).*(longSignal(25:1200*samplesPerCode+24)));
longSignal = double(longSignal(1:20*samplesPerCode));
%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq   = zeros(1, 37);
% C/A code phases of detected signals
acqResults.codePhase  = zeros(1, 37);
% Correlation peak ratios of the detected signals
acqResults.peakMetric = zeros(1, 37);
acqResults.codeDopple = zeros(1, 37);
caCodes = single(makeCaTable(settings,10,0));
msOfTimeIndex = unique(round(102.3./(-3:0.01:3)),'stable');
msOfTimeIndex((abs(msOfTimeIndex)>1000)&(~isinf(msOfTimeIndex)))=[];
dopplerCode = (102.3)./msOfTimeIndex;
result  =  zeros(37,length(msOfTimeIndex),samplesPerCode);
fprintf('(');
for i = 1:length(msOfTimeIndex)
    if isinf(msOfTimeIndex(i))
        corrSignal = sum(reshape(signal(1:1000*samplesPerCode),samplesPerCode,1000),2)';
    elseif msOfTimeIndex(i) < 0
        m = - msOfTimeIndex(i);
        corrTime    = abs(floor(1000/(m)));
        lenOfTime   = m * samplesPerCode + 1;
        corrSignal  = sum(reshape(signal(1:corrTime*lenOfTime),lenOfTime,corrTime),2);
        corrSignal  = sum(reshape(corrSignal(1:(end-1)),samplesPerCode,m),2)';
    else
        corrTime    = floor(1000/(msOfTimeIndex(i)));
        lenOfTime   = msOfTimeIndex(i) * samplesPerCode - 1;
        corrSignal  = sum(reshape(signal(1:corrTime*lenOfTime),lenOfTime,corrTime),2);
        corrSignal  = sum(reshape([corrSignal;0],samplesPerCode,msOfTimeIndex(i)),2)';
    end
    signalFreqDom = fft(single(corrSignal));
    result(:,i,:) = bsxfun(@times,signalFreqDom,fftCode);
end
result = abs(ifft(result,samplesPerCode,3));
        %% Correlate signals ======================================================
for PRN = settings.acqSatelliteList      
    %% Look for correlation peaks in the results ==============================
    [~ ,codeDoppleIndex] = max(max(result(PRN,:,:), [], 3));
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
    %--- Find code phase of the same correlation peak ---------------------
    [peakSize ,codePhase ]  = max(max(result(PRN,:,:))); 
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip = round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;
    
    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 1
        codePhaseRange = excludeRangeIndex2 : ...
            (samplesPerCode + excludeRangeIndex1);
    elseif excludeRangeIndex2 > samplesPerCode
        codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
            excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
            excludeRangeIndex2 : samplesPerCode];
    end
    
    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(result(PRN,codeDoppleIndex, codePhaseRange));
    
    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(PRN) = peakSize/secondPeakSize;
    acqResults.codePhase (PRN) = codePhase;
    acqResults.codeDopple(PRN) = dopplerCode(codeDoppleIndex);
    % If the result is above threshold, then there is a signal ...
    if  acqResults.peakMetric(PRN) > 1.5
        %% Fine resolution frequency search =======================================
        %--- Indicate PRN number of the detected signal -------------------
        fprintf('%02d ', PRN);
        xCarrier = ...
            longSignal(codePhase:codePhase + 10*samplesPerCode - 1) ...
            .* caCodes(PRN,:);
        xCarrier = reshape(xCarrier,samplesPerCode,10);
        %--- Find the next highest power of two and increase by 8x --------
        cztxc=abs(czt(xCarrier,1526,exp(-1i*2*pi/settings.samplingFreq),...
            exp(1i*2*pi*((acqResults.codeDopple(PRN)-1)*763+settings.IF)/settings.samplingFreq)));
        cztxc=sum(cztxc,2);
        cztFreqBins=(acqResults.codeDopple(PRN)-1)*763+settings.IF:1:(acqResults.codeDopple(PRN)+1)*763+settings.IF;
        [~, cztMaxIndex] = max(cztxc);
        %--- Save properties of the detected satellite signal -------------
        acqResults.carrFreq(PRN)  = cztFreqBins(cztMaxIndex);      
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end
end% for PRN = satelliteList
%=== Acquisition is over ==================================================
fprintf(')\n');
