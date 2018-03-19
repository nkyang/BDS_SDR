function [p_fa,p_d] = acquisition1_monteCarlo(data,fftCode, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 11 ms of raw signal from the front-end
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Darius Plausinaitis and Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on Peter Rinder and Nicolaj Bertelsen
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $

%% Initialization =========================================================

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
samplesPerCodeChip = round(settings.samplingFreq / settings.codeFreqBasis);
% Generate all C/A codes and sample them according to the sampling freq.
% caCodesTable = makeCaTable(settings,1,0);
% NHcode = [-1,-1,-1,-1,-1,1,-1,-1,1,1,-1,1,-1,1,-1,-1,1,1,1,-1];
% longNHcacode = kron(NHcode,caCodesTable);
%--- Initialize arrays to speed up the code -------------------------------
peakMetric = zeros(37, 1);
result     = zeros(37,33, samplesPerCode);
trueSat    = [1,2,3,4,5,6,7,9,13,14,17,32];
wrongSat   = setdiff(1:37,trueSat);
thresHold  = 1.5:0.01:2.5;
nFalse     = zeros(1,length(thresHold));
nCorre     = zeros(37,length(thresHold));
for num = 1:10000
    rStart = randi(30000*samplesPerCode,1,1);
    longSignal=data(rStart:(1100 * samplesPerCode+rStart-1));
    % Create two 1msec vectors of data to correlate with and one with zero DC
    signal = longSignal(1 : 1030 * samplesPerCode).*longSignal(25:(1030 * samplesPerCode + 24));
    for dopplerCodeIndex = 1:33
        dopplerCode = (dopplerCodeIndex - 17)*0.5;
        if dopplerCode == 0
            corrSignal = sum(reshape(signal(1:1000*samplesPerCode),samplesPerCode,1000),2)';
        elseif dopplerCode > 0
            msOfTime    = round(102.3/dopplerCode);
            corrTime    = round(1000/(msOfTime));
            lenOfTime   = msOfTime * samplesPerCode - 1;
            corrSignal  = sum(reshape(signal(1:corrTime*lenOfTime),lenOfTime,corrTime),2);
            corrSignal  = sum(reshape([corrSignal;0],samplesPerCode,msOfTime),2)';
        else
            msOfTime    = abs(round(102.3/dopplerCode));
            corrTime    = round(1000/(msOfTime));
            lenOfTime   = msOfTime * samplesPerCode + 1;
            corrSignal  = sum(reshape(signal(1:corrTime*lenOfTime),lenOfTime,corrTime),2);
            corrSignal  = sum(reshape(corrSignal(1:(end-1)),samplesPerCode,msOfTime),2)';
        end
        signalFreqDom = fft(corrSignal);
        result(:,dopplerCodeIndex,:) = bsxfun(@times,signalFreqDom,fftCode);
    end
    result = abs(ifft(result, samplesPerCode,3)).^2;
    for PRN = settings.acqSatelliteList
        %% Look for correlation peaks in the results ==============================
        [~ ,codeDoppleIndex] = max(max(result(PRN,:,:), [], 3));
        % Find the highest peak and compare it to the second highest peak
        % The second peak is chosen not closer than 1 chip to the highest peak
        %--- Find code phase of the same correlation peak ---------------------
        [peakSize ,codePhase ]  = max(max(result(PRN,:,:)));
        %--- Find 1 chip wide C/A code phase exclude range around the peak ----
        
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
        peakMetric(PRN) = peakSize/secondPeakSize;
        
    end
    %=== Acquisition is over ==================================================
    % fprintf(')\n');
    a = bsxfun(@ge,peakMetric,thresHold);
    nCorre = nCorre + a;
    nFalse = nFalse + sum(a(wrongSat,:));
end
p_fa = nFalse /(num  * (37 - length(trueSat)));
p_d  = nCorre / num;
save ROC_1s.mat p_fa p_d

