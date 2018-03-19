function [eph, SOW] = ephemeris_D1(bits)
%Function decodes ephemerides and TOW from the given bit stream. The stream
%(array) in the parameter BITS must contain 1500 bits. The first element in
%the array must be the first bit of a subframe. The subframe ID of the
%first subframe in the array is not important.
%
%Function does not check parity!
%
%[eph, TOW] = ephemeris(bits, D30Star)
%
%   Inputs:
%       bits        - bits of the navigation messages (5 subframes).
%                   Type is character array and it must contain only
%                   characters '0' or '1'.
%       D30Star     - The last bit of the previous nav-word. Refer to the
%                   GPS interface control document ICD (IS-GPS-200D) for
%                   more details on the parity checking. Parameter type is
%                   char. It must contain only characters '0' or '1'.
%   Outputs:
%       TOW         - Time Of Week (TOW) of the first sub-frame in the bit
%                   stream (in seconds)
%       eph         - SV ephemeris

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Darius Plausinaitis and Kristin Larson
% Written by Darius Plausinaitis and Kristin Larson
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
%$Id: ephemeris.m,v 1.1.2.7 2006/08/14 11:38:22 dpl Exp $


%% Check if there is enough data ==========================================
% if length(bits) < 1500
%     error('The parameter BITS must contain 1500 bits!');
% end

%% Check if the parameters are strings ====================================


% Pi used in the GPS coordinate system
bdsPi = 3.1415926535898;

preamble_bits = [1,1,1,0,0,0,1,0,0,1,0];
eph=struct('health',[],'AODC',[],'URAI',[],'weekNumber',[],'t_oc',[],'T_GD',[],'alpha0',[],'alpha1',[],...
    'alpha2',[],'alpha3',[],'beta0',[],'beta1',[],'beta2',[],'beta3',[],'a0',[],'a1',[],...
    'a2',[],'AODE',[],'deltan',[],'C_uc',[],'M_0',[],'C_us',[],'e',[],'sqrtA',[],...
    'C_ic',[],'C_is',[],'t_oe',[],'i_0',[],'C_rc',[],'C_rs',[],'omega_dot',[],'omega_0',[],...
    'omega',[],'IDOT',[]);
%% Decode all 5 sub-frames ================================================
for i = 1:1:5
    
    %--- "Cut" one sub-frame's bits ---------------------------------------
    subframe = bits(300*(i-1)+1 : 300*i);
    if ((preamble_bits*subframe(1:11)) < 5)
        subframe = ~subframe;
    end
    [subframe(16:30),~] = BCH(subframe(16:30)');
    %--- Correct polarity of the data bits in all 10 words ----------------
    for j = 2:10
        subframe(30*(j-1)+1 : 30*j) = deinter(subframe(30*(j-1)+1 : 30*j));
    end
    subframe = num2str(subframe+0)';
    %--- Decode the sub-frame id ------------------------------------------
    % For more details on sub-frame contents please refer to GPS IS.
    subframeID = bin2dec(subframe(16:18));
    %--- Decode sub-frame based on the sub-frames id ----------------------
    % The task is to select the necessary bits and convert them to decimal
    % numbers. For more details on sub-frame contents please refer to GPS
    % ICD (IS-GPS-200D).
    switch subframeID 
        % It contains WN, SV clock corrections, health and accuracy
        case 1  %--- It is subframe 1 -------------------------------------
            SOW             = bin2dec([subframe(19:26) subframe(31:42)]);
            eph.health      = bin2dec(subframe(43));
            eph.AODC        = bin2dec(subframe(44:48));
            eph.URAI        = bin2dec(subframe(49:52));
            eph.weekNumber  = bin2dec(subframe(61:73));
            eph.t_oc        = bin2dec([subframe(74:82) subframe(91:98)]) * 2^3;
            eph.T_GD        = twosComp2dec(subframe(99:109)) * 10^(-10);
            eph.alpha0      = twosComp2dec(subframe(127:134)) * 2^(-30);
            eph.alpha1      = twosComp2dec(subframe(135:142)) * 2^(-27);
            eph.alpha2      = twosComp2dec(subframe(151:158)) * 2^(-24);
            eph.alpha3      = twosComp2dec(subframe(159:166)) * 2^(-24);
            eph.beta0       = twosComp2dec([subframe(167:172) subframe(181:182)]) * 2^11;
            eph.beta1       = twosComp2dec(subframe(183:190)) * 2^14;
            eph.beta2       = twosComp2dec(subframe(191:198)) * 2^16;
            eph.beta3       = twosComp2dec([subframe(199:202) subframe(211:214)]) * 2^16;
            eph.a2          = twosComp2dec(subframe(215:225)) * 2^(-66);
            eph.a0          = twosComp2dec([subframe(226:232) subframe(241:257)]) * 2^(-33);
            eph.a1          = twosComp2dec([subframe(258:262) subframe(271:287)]) * 2^(-50);
            eph.AODE        = bin2dec(subframe(288:292));
        case 2  %--- It is subframe 2 -------------------------------------
            eph.deltan      = twosComp2dec([subframe(43:52) subframe(61:66)]) * 2^(-43) * bdsPi;
            eph.C_uc        = twosComp2dec([subframe(67:82) subframe(91:92)]) * 2^(-31);
            eph.M_0         = twosComp2dec([subframe(93:112) subframe(121:132)]) * 2^(-31) * bdsPi;
            eph.e           = twosComp2dec([subframe(133:142) subframe(151:172)]) * 2^(-33);
            eph.C_us        = twosComp2dec(subframe(181:198) ) * 2^(-31);
            eph.C_rc        = twosComp2dec([subframe(199:202) subframe(211:224)]) * 2^(-6);
            eph.C_rs        = twosComp2dec([subframe(225:232) subframe(241:250)]) * 2^(-6);
            eph.sqrtA       = bin2dec([subframe(251:262) subframe(271:290)]) * 2^(-19);
            t_oe_MSB        = subframe(291:292);
        case 3  %--- It is subframe 3 -------------------------------------
            eph.t_oe        = bin2dec([t_oe_MSB subframe(43:52) subframe(61:65)]) * 2^3;
            eph.i_0         = twosComp2dec([subframe(66:82) subframe(91:105)]) * 2^(-31) * bdsPi;
            eph.C_ic        = twosComp2dec([subframe(106:112) subframe(121:131)]) * 2^(-31);
            eph.omega_dot   = twosComp2dec([subframe(132:142) subframe(151:163)]) * 2^(-43) * bdsPi;
            eph.C_is        = twosComp2dec([subframe(164:172) subframe(181:189)]) * 2^(-31);
            eph.IDOT        = twosComp2dec([subframe(190:202) subframe(211)]) * 2^(-43) * bdsPi;
            eph.omega_0     = twosComp2dec([subframe(212:232) subframe(241:251)]) * 2^(-31) * bdsPi;
            eph.omega       = twosComp2dec([subframe(252:262) subframe(271:291)]) * 2^(-31) * bdsPi;
%         case 4
%             eph.pnum4=bin2dec(subframe(44:50));
%         case 5
%             eph.pnum5=bin2dec(subframe(44:50));
    end % switch subframeID ...
    
end % for all 5 sub-frames ...



%% Compute the time of week (TOW) of the first sub-frames in the array ====
% Also correct the TOW. The transmitted TOW is actual TOW of the next
% subframe and we need the TOW of the first subframe in this data block
% (the variable subframe at this point contains bits of the last subframe).

