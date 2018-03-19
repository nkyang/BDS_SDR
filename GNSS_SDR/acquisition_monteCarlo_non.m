function [p_fa,p_d] = acquisition_monteCarlo_diff(longSignal, settings)
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
% Find sampling period
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave
% phasePoints_D1 = (0 : (20*samplesPerCode-1)) * 2 * pi * ts;
phasePoints_D2 = (0 : (2*samplesPerCode-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band (500Hz steps)
numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;

% Generate all C/A codes and sample them according to the sampling freq.
caCodesTable = makeCaTable(settings,2,0);
caCodeFreqDom = reshape(conj(fft(caCodesTable,2*samplesPerCode,2)),[37,1,2*samplesPerCode]);
%--- Initialize arrays to speed up the code -------------------------------
% Carrier frequencies of the frequency bins
frqBins   = settings.IF - (settings.acqSearchBand/2) * 1000 + 500 * (0:numberOfFrqBins - 1);
sinCarr   = sin(bsxfun(@times,frqBins',phasePoints_D2));
cosCarr   = cos(bsxfun(@times,frqBins',phasePoints_D2));

trueSat    = [1,2,3,4,5,6,7,9,13,14,17,32];
wrongSat   = setdiff(1:37,trueSat);
thresHold  = 1.2:0.001:1.7;
nFalse     = zeros(1,length(thresHold));
nCorre     = zeros(37,length(thresHold));
peakMetric   = zeros(37,1);
for num = 1:10000
    
    rStart = randi(3e7,1,1);
    % Create two 1msec vectors of data to correlate with and one with zero DC
    signal1 = double(longSignal(rStart : 2 * samplesPerCode + rStart - 1));
    signal2 = double(longSignal(2 * samplesPerCode + rStart : 4*samplesPerCode + rStart - 1));
    %--- "Remove carrier" from the signal -----------------------------
    I1      = bsxfun(@times,sinCarr , signal1);
    Q1      = bsxfun(@times,cosCarr , signal1);
    I2      = bsxfun(@times,sinCarr , signal2);
    Q2      = bsxfun(@times,cosCarr , signal2);
    %--- Multiplication in the frequency domain (correlation in time
    %domain)
    IQ1 = fft(I1 + 1j*Q1,2*samplesPerCode,2);
    IQ2 = fft(I2 + 1j*Q2,2*samplesPerCode,2);
    IQ1 = reshape(IQ1,[1,numberOfFrqBins,2*samplesPerCode]);
    IQ2 = reshape(IQ2,[1,numberOfFrqBins,2*samplesPerCode]);
    convCodeIQ1 = bsxfun(@times,IQ1,caCodeFreqDom);
    convCodeIQ2 = bsxfun(@times,IQ2,caCodeFreqDom);
    %--- Perform inverse DFT and store correlation results ------------
    acqRes1 = ifft(convCodeIQ1,2*samplesPerCode,3) ;
    acqRes2 = ifft(convCodeIQ2,2*samplesPerCode,3) ;
    results = abs(acqRes1)+abs(acqRes2);
    %% Look for correlation peaks in the results ===========[===================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
    for PRN = settings.acqSatelliteList
        %% Look for correlation peaks in the results ==============================
        [~ ,codeDoppleIndex] = max(max(results(PRN,:,:), [], 3));
        % Find the highest peak and compare it to the second highest peak
        % The second peak is chosen not closer than 1 chip to the highest peak
        %--- Find code phase of the same correlation peak ---------------------
        [peakSize ,codePhase ]  = max(max(results(PRN,:,:)));
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
        secondPeakSize = max(results(PRN,codeDoppleIndex, codePhaseRange));
        %--- Store result -----------------------------------------------------
        peakMetric(PRN) = peakSize/secondPeakSize;
    end
    a = bsxfun(@ge,peakMetric,thresHold);
    nCorre = nCorre + a;
    nFalse = nFalse + sum(a(wrongSat,:));
end
p_fa = nFalse /(num  * (37 - length(trueSat)));
p_d  = nCorre / num;
save ROC_2+2non.mat p_fa p_d

