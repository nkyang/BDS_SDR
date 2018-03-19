function acqResults = acquisition(longSignal, settings)
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

% Find sampling period
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave
% phasePoints_D1 = (0 : (20*samplesPerCode-1)) * 2 * pi * ts;
phasePoints_D2 = (0 : (2*samplesPerCode-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band (500Hz steps)
numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;

% Generate all C/A codes and sample them according to the sampling freq.
caCodesTable = single(makeCaTable(settings,1,0));
% NHcode = [-1,-1,-1,-1,-1,1,-1,-1,1,1,-1,1,-1,1,-1,-1,1,1,1,-1];
% longNHcacode = kron(NHcode,caCodesTable);
%--- Initialize arrays to speed up the code -------------------------------
% Search results of all frequency bins and code shifts (for one satellite)
% results_D1     = zeros(numberOfFrqBins, 20*samplesPerCode);
results_D2  = zeros(numberOfFrqBins, 2*samplesPerCode);
% Carrier frequencies of the frequency bins
frqBins   = zeros(1, numberOfFrqBins);


%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 37);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 37);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 37);

fprintf('(');

% Create two 1msec vectors of data to correlate with and one with zero DC
% signal20 = longSignal(1 : 20*samplesPerCode);
signal1 = longSignal(1 : 2 * samplesPerCode);
signal2 = longSignal(2 * samplesPerCode + 1 : 4 * samplesPerCode );
signalDC = longSignal - mean(longSignal);
% Perform search for all listed PRN numbers ...
for PRN = settings.acqSatelliteList
        caCodeFreqDom1 = conj(fft([caCodesTable(PRN, :),caCodesTable(PRN, :)]));
%         caCodeFreqDom2 = conj(fft([caCodesTable(PRN, 9996:10000),caCodesTable(PRN, :),caCodesTable(PRN, 1:9995)]));
        
        for frqBinIndex = 1:numberOfFrqBins
            
            %--- Generate carrier wave frequency grid (0.5kHz step) -----------
            frqBins(frqBinIndex) = settings.IF - ...
                (settings.acqSearchBand/2) * 400 + ...
                0.2e3 * (frqBinIndex - 1);
            
            %--- Generate local sine and cosine -------------------------------
            sinCarr = sin(frqBins(frqBinIndex) * phasePoints_D2);
            cosCarr = cos(frqBins(frqBinIndex) * phasePoints_D2);
            
            %--- "Remove carrier" from the signal -----------------------------
            I1      = sinCarr .* signal1;
            Q1      = cosCarr .* signal1;
            I2      = sinCarr .* signal2;
            Q2      = cosCarr .* signal2;  
            %--- Multiplication in the frequency domain (correlation in time
            %domain)
            convCodeIQ1 = fft(I1 + 1j*Q1) .* caCodeFreqDom1; 
            convCodeIQ2 = fft(I2 + 1j*Q2) .* caCodeFreqDom1;
            %--- Perform inverse DFT and store correlation results ------------
            acqRes1 = ifft(convCodeIQ1) ;
            acqRes2 = ifft(convCodeIQ2) ;
            results_D2(frqBinIndex, :) = abs(acqRes1.*acqRes2);
        end % frqBinIndex = 1:numberOfFrqBins
        
        %% Look for correlation peaks in the results ==============================
        % Find the highest peak and compare it to the second highest peak
        % The second peak is chosen not closer than 1 chip to the highest peak
        
        %--- Find the correlation peak and the carrier frequency --------------
        [~ ,frequencyBinIndex] = max(max(results_D2, [], 2));
        
        %--- Find code phase of the same correlation peak ---------------------
        [peakSize ,codePhase] = max(max(results_D2));
         %--- Find 1 chip wide C/A code phase exclude range around the peak ----
        samplesPerCodeChip = round(settings.samplingFreq / settings.codeFreqBasis);
        excludeRangeIndex1 = codePhase - samplesPerCodeChip;
        excludeRangeIndex2 = codePhase + samplesPerCodeChip;
        
        %--- Correct C/A code phase exclude range if the range includes array
        %boundaries
        if excludeRangeIndex1 < 1
            codePhaseRange = excludeRangeIndex2 : (samplesPerCode + excludeRangeIndex1);         
        elseif excludeRangeIndex2 > samplesPerCode
            codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : excludeRangeIndex1;
        else
            codePhaseRange = [1:excludeRangeIndex1, excludeRangeIndex2 : samplesPerCode];
        end
        
        %--- Find the second highest correlation peak in the same freq. bin ---
        secondPeakSize = max(results_D2(frequencyBinIndex, codePhaseRange));
             
        %--- Store result -----------------------------------------------------
        acqResults.peakMetric(PRN) = (peakSize/secondPeakSize);
        % If the result is above threshold, then there is a signal ...   
        if (peakSize/secondPeakSize) > 2.2
%         if (peakSize/mean(results_D2(frequencyBinIndex,:))) > 12.5
            %% Fine resolution frequency search =======================================
            
            %--- Indicate PRN number of the detected signal -------------------
            fprintf('%02d ', PRN);

            %--- Remove C/A code modulation from the original signal ----------
            % (Using detected C/A code phase)
            xCarrier = ...
                signalDC(codePhase:(codePhase + 10*samplesPerCode-1)) .* kron(ones(1,10),caCodesTable(PRN, :));
            xCarrier = reshape(xCarrier,samplesPerCode,10);
            cztxc=abs(czt(xCarrier,401,exp(-1i*2*pi/settings.samplingFreq),exp(1i*2*pi*frqBins(frequencyBinIndex-1)/settings.samplingFreq)));
            cztxc=sum(cztxc,2);
            cztFreqBins=frqBins(frequencyBinIndex-1):1:frqBins(frequencyBinIndex+1);
%             cztxc=abs(czt(xCarrier,20001,exp(-1i*2*pi*20000/(20001*1e7)),exp(1i*2*pi*0.249)));
%             cztFreqBins=(settings.IF-1e4):1:(settings.IF+1e4);
            [~, cztMaxIndex] = max(cztxc);
            acqResults.carrFreq(PRN)  = cztFreqBins(cztMaxIndex);
            acqResults.codePhase(PRN) = codePhase;
        else
            %--- No signal with this PRN --------------------------------------
            fprintf('. ');
        end   % if (peakSize/secondPeakSize) > settings.acqThreshold
end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
