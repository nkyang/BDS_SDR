function acqResults = acquisition2(rawSignal,cumSignal,fftCode,settings, interTime)
%% Initialization =========================================================
% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
%--- Find 1 chip wide C/A code phase exclude range around the peak ----
samplesPerCodeChip = round(settings.samplingFreq / settings.codeFreqBasis);
% Create two 1msec vectors of data to correlate with and one with zero DC
rawSignal = double(rawSignal(1:101*samplesPerCode));
satNum = length(settings.acqSatelliteList);
%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq   = zeros(1, satNum);
% ranging code phases of detected signals
acqResults.codePhase  = zeros(1, satNum);
% Dopple shift of the ranging code
acqResults.codeDopple = zeros(1, satNum);
% Correlation peak ratios of the detected signals
acqResults.peakMetric = zeros(1, satNum);
% intergral time
% interTime = 12000;
% Correlatoin result of the detected signals
codeDopple = -4:0.005:4;
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
blockSize  = round(settings.codeFreqBasis ./ codeDopple);
blockSize((blockSize>interTime * samplesPerCode)&(~isinf(blockSize))) = [];
codeDopple = settings.codeFreqBasis ./ blockSize;
% blockSize  = abs(blockSize);
signalFreqDom = fft(cumSignal,samplesPerCode,2);
clear cumSignal
fprintf('(');
    %% Correlate signals ======================================================
for ii = 1:satNum
    a = bsxfun(@times, signalFreqDom, fftCode(ii,:));
    result = abs(ifft(a, samplesPerCode, 2)).^2;
    [~ ,codeDoppleIndex ] = max(max(result,[],2),[],1);
    codeFreqShift         = codeDopple(codeDoppleIndex);
    [peakSize ,codePhase] = max(max(result,[],1),[],2);
    id1 = codePhase - samplesPerCodeChip;
    id2 = codePhase + samplesPerCodeChip;
    A   = 1:samplesPerCode;
    codePhaseRange = A((A<id1)&(A>id2-samplesPerCode)|((A>id2)&(A<id1+samplesPerCode)));
    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(result(codeDoppleIndex,codePhaseRange));
    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(ii) = peakSize/secondPeakSize;
    if  acqResults.peakMetric(ii) > 2.2
        %% Fine resolution frequency search =======================================
        %--- Indicate PRN number of the detected signal -------------------
        PRN = settings.acqSatelliteList(ii);
        fprintf('%02d ', PRN);
        caCodes  = makeCaTable(settings,100,codeFreqShift,PRN);
        xCarrier = rawSignal(codePhase:codePhase + length(caCodes) - 1) .* caCodes;
%         xCarrier = reshape(xCarrier,samplesPerCode,20);
        cztxc    = abs(czt(xCarrier,2001,exp(-1i*2*pi*1/settings.samplingFreq),...
            exp(1i*2*pi*(codeFreqShift*763-1000+settings.IF)/settings.samplingFreq)));
%         cztxc    = sum(cztxc,2);
        cztFreqBins = (codeFreqShift*763-1000+settings.IF):1:...
            (codeFreqShift*763+1000+settings.IF);
        [~, cztMaxIndex] = max(cztxc);
        figure();
        plot(cztFreqBins,cztxc);
        %--- Save properties of the detected satellite signal -------------
        acqResults.carrFreq(ii)  = cztFreqBins(cztMaxIndex);      
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end
    acqResults.codePhase(ii)  = codePhase;
    acqResults.codeDopple(ii) = codeDopple(codeDoppleIndex);
    acqResults.result(ii,:,:) = result;
end    
% %% Look for correlation peaks in the results ==============================
% [~ ,codeDoppleIndex ] = max(max(result,[],3),[],2);
% [peakSize ,codePhase] = max(max(result,[],2),[],3); 
% acqResults.codePhase  = codePhase;
% acqResults.codeDopple = codeDopple(codeDoppleIndex);
% acqResults.result     = result;
% for ii = 1:satNum 
%     id1 = codePhase(ii) - samplesPerCodeChip;
%     id2 = codePhase(ii) + samplesPerCodeChip;
%     A   = 1:samplesPerCode;
%     codePhaseRange = A((A<id1)&(A>id2-samplesPerCode)|((A>id2)&(A<id1+samplesPerCode)));
%     %--- Find the second highest correlation peak in the same freq. bin ---
%     secondPeakSize = max(result(ii,codeDoppleIndex(ii),codePhaseRange));    
%     %--- Store result -----------------------------------------------------
%     acqResults.peakMetric(ii) = peakSize/secondPeakSize;
%     % If the result is above threshold, then there is a signal ...
%     if  acqResults.peakMetric(ii) > 2.2
%         %% Fine resolution frequency search =======================================
%         %--- Indicate PRN number of the detected signal -------------------
%         PRN = settings.acqSatelliteList(ii);
%         fprintf('%02d ', PRN);
%         caCodes  = makeCaTable(settings,30,acqResults.codeDopple(ii),PRN);
%         xCarrier = rawSignal(codePhase(ii):codePhase(ii) + 30*samplesPerCode - 1) .* caCodes;
%         xCarrier = reshape(xCarrier,samplesPerCode,30);
%         cztxc    = abs(czt(xCarrier,201,exp(-1i*2*pi*10/settings.samplingFreq),...
%             exp(1i*2*pi*(acqResults.codeDopple(ii)*763-1000+settings.IF)/settings.samplingFreq)));
%         cztxc    = sum(cztxc,2);
%         cztFreqBins = (acqResults.codeDopple(ii)*763-1000+settings.IF):10:...
%             (acqResults.codeDopple(ii)*763+1000+settings.IF);
%         [~, cztMaxIndex] = max(cztxc);
%         figure();
%         plot(cztFreqBins,cztxc);
%         %--- Save properties of the detected satellite signal -------------
%         acqResults.carrFreq(ii)  = cztFreqBins(cztMaxIndex);      
%     else
%         %--- No signal with this PRN --------------------------------------
%         fprintf('. ');
%     end
% end% for PRN = satelliteList
%=== Acquisition is over ==================================================
fprintf(')\n');
