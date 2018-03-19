function [eph, SOW] = ephemeris_D2(bits)
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
for i = 1:1:50
    
    %--- "Cut" one sub-frame's bits ---------------------------------------
    subframe = bits(300*(i-1)+1 : 300*i);
    if ((preamble_bits*subframe(1:11)) < 5 )
        subframe = ~subframe;
    end
    %--- Correct polarity of the data bits in all 10 words ----------------
    [subframe(16:30),~] = BCH(subframe(16:30)');
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
    SOW             = bin2dec([subframe(19:26) subframe(31:42)]) ;
    if subframeID==1 %--- It is subframe 1 -------------------------------------
        % It contains WN, SV clock corrections, health and accuracy
        
        pnum            = bin2dec(subframe(43:46));
        switch pnum
            case 1
                eph.health      = bin2dec(subframe(47));
                eph.AODC        = bin2dec(subframe(48:52));
                eph.URAI        = bin2dec(subframe(61:64));
                eph.weekNumber  = bin2dec(subframe(65:77));
                eph.t_oc        = bin2dec([subframe(78:82) subframe(91:102)]) * 2^3;
                eph.T_GD        = twosComp2dec(subframe(103:112)) * 10^(-10);
                SOW             = bin2dec([subframe(19:26) subframe(31:42)]) ;
            case 2
                eph.alpha0      = twosComp2dec([subframe(47:52) subframe(61:62)]) * 2^(-30);
                eph.alpha1      = twosComp2dec(subframe(63:70)) * 2^(-27);
                eph.alpha2      = twosComp2dec(subframe(71:78)) * 2^(-24);
                eph.alpha3      = twosComp2dec([subframe(79:82) subframe(91:94)]) * 2^(-24);
                eph.beta0       = twosComp2dec(subframe(95:102)) * 2^11;
                eph.beta1       = twosComp2dec(subframe(103:110)) * 2^14;
                eph.beta2       = twosComp2dec([subframe(111:112) subframe(121:126)]) * 2^16;
                eph.beta3       = twosComp2dec(subframe(127:134)) * 2^16;
                SOW             = bin2dec([subframe(19:26) subframe(31:42)]) ;
            case 3
                eph.a0          = twosComp2dec([subframe(101:112) subframe(121:132)]) * 2^(-33);
                a1_MSB          = subframe(133:136);
                SOW             = bin2dec([subframe(19:26) subframe(31:42)]) ;
            case 4
                eph.a1          = twosComp2dec([a1_MSB subframe(47:52) subframe(61:72)]) * 2^(-50);
                eph.a2          = twosComp2dec([subframe(73:82) subframe(91)]) * 2^(-66);
                eph.AODE        = bin2dec(subframe(92:96));
                eph.deltan      = twosComp2dec(subframe(97:112)) * 2^(-43) * bdsPi;
                C_uc_MSB        = subframe(121:134);
                SOW             = bin2dec([subframe(19:26) subframe(31:42)]) ;
            case 5
                eph.C_uc        = twosComp2dec([C_uc_MSB subframe(47:50)]) * 2^(-31);
                eph.M_0         = twosComp2dec([subframe(51:52) subframe(61:82) subframe(91:98)]) * 2^(-31) * bdsPi;
                eph.C_us        = twosComp2dec([subframe(99:112) subframe(121:124)]) * 2^(-31);
                e_MSB           = subframe(125:134);
                SOW             = bin2dec([subframe(19:26) subframe(31:42)]) ;
            case 6
                eph.e           = bin2dec([e_MSB subframe(47:52) subframe(61:76)]) * 2^(-33);
                eph.sqrtA       = bin2dec([subframe(77:82) subframe(91:112) subframe(121:124)]) * 2^(-19);
                C_ic_MSB        = subframe(125:134);
                SOW             = bin2dec([subframe(19:26) subframe(31:42)]) ;
            case 7
                eph.C_ic        = twosComp2dec([C_ic_MSB subframe(47:62) subframe(61:62)]) * 2^(-31);
                eph.C_is        = twosComp2dec(subframe(63:80)) * 2^(-31);
                eph.t_oe        = bin2dec([subframe(81:82) subframe(91:105)]) * 2^3;
                i_0_MSB         = [subframe(106:112) subframe(121:134)];
                SOW             = bin2dec([subframe(19:26) subframe(31:42)]) ;
            case 8
                eph.i_0         = twosComp2dec([i_0_MSB subframe(47:52) subframe(61:65)]) * 2^(-31) * bdsPi;
                eph.C_rc        = twosComp2dec([subframe(66:82) subframe(91)]) * 2^(-6);
                eph.C_rs        = twosComp2dec(subframe(92:109)) * 2^(-6);
                omega_dot_MSB   = [subframe(110:112) subframe(121:136)];
                SOW             = bin2dec([subframe(19:26) subframe(31:42)]) ;
            case 9
                eph.omega_dot   = twosComp2dec([omega_dot_MSB subframe(47:51)]) * 2^(-43) * bdsPi;
                eph.omega_0     = twosComp2dec([subframe(52) subframe(61:82) subframe(91:99)]) * 2^(-31) * bdsPi;
                omega_MSB       = [subframe(100:112) subframe(121:134)];
                SOW             = bin2dec([subframe(19:26) subframe(31:42)]) ;
            case 10
                eph.omega       = twosComp2dec([omega_MSB subframe(47:51)]) * 2^(-31) * bdsPi;
                eph.IDOT        = twosComp2dec([subframe(52) subframe(61:73)]) * 2^(-43) * bdsPi;
                SOW             = bin2dec([subframe(19:26) subframe(31:42)]) ;
        end %switch pnum
        
    end % switch subframeID ...
    
end % for all 5 sub-frames ...
%                   eph.a1          = twosComp2dec([eph.a1_MSB e.pha1_LSB]) * 2^(-50);
%                   eph.C_uc        = twosComp2dec([eph.C_uc_MSB eph.C_uc_LSB]) * 2^(-31);
%                   eph.e           = bin2dec([eph.e_MSB eph.e_LSB]) * 2^(-33);
%                   eph.C_ic        = twosComp2dec([eph.C_ic_MSB eph.C_ic_LSB]) * 2^(-31);
%                   eph.i_0         = twosComp2dec([eph.i_0_MSB eph.i_0_LSB]) * 2^(-31) * bdsPi;
%                   eph.omega_dot   = twosComp2dec([eph.omega_dot_MSB eph.omega_dot_LSB]) * 2^(-43) * bdsPi;
%                   eph.omega       = twosComp2dec([eph.omega_MSB eph.omega_LSB]) * 2^(-31) * bdsPi;


%% Compute the time of week (TOW) of the first sub-frames in the array ====
% Also correct the TOW. The transmitted TOW is actual TOW of the next
% subframe and we need the TOW of the first subframe in this data block
% (the variable subframe at this point contains bits of the last subframe).

