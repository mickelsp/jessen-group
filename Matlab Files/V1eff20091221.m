function V120091221=V1eff20091221(Podt,temp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% V1eff20091221 = V1eff20091221(PowerODT, temp) where Powers are in Watts and temp in
%% Kelvin based on the waists extracted from the trap frequencies measured
%% on Dec 17, 2009 and waists are redetermined on Dec 21, 2009, referring to Page 46 of Mi's notebook.
%% V1eff is defined as the effective volume V1 (refer to Eq.3 in our ODT model paper arXiv:0905.2223)
%% We have n_peak * V1 = N, where N is the total number of atoms and n_peak
%% is the peak density of atoms in ODT.
%% We built lookup table for V1eff as the function of the ODT
%% power and equilibrium temperature using Mathematica 6.0 in advance.
%% Then create the interpolation function for V1eff here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Note1: the range of variables in lookup table are 
%% PowerODT : 1.078 - 9.8 W
%% equlibrium Temp : 0.00164 - 32.6 uK

%% Although the Cubic spline interpolation ('spline') used here
%% performs extrapolation for out of range values, I recommend to input
%% variable within the lookup table range mentioned above. If it is
%% necessary
%% to get relatively reasonable values for V1eff out of the range, you'd 
%% better build more lookup table within the range you may be interested in for the
%% interpolation.


directory=char('/Users/work/Documents/Analysis/M files/V1ODTLookupTable12212009/');
    

if ((Podt > 1.47) & (Podt <= 9.8)) %Interperation from lookup table 1 where the power valid range is 1.47W to 9.8W
   SampNumPowerODT=9; % There are 9 sample ODT power values: 1.47,  ... , 9.8 W
   SampNumTemp=200; % There are 100 sample equlibrium temperature values: 1*32.6uK/200, 2*32.6uK/200, ... , 100*32.6uK/200
   sumdatafile=strcat(directory,'Table1V1effDec212009.dat');
   [PowerODTTemp TempODTDimpleTemp VeffDataTemp]=textread(sumdatafile, '%f%f%f','commentstyle','matlab');
   Veff=zeros(SampNumPowerODT,SampNumTemp);
            for i=1:SampNumPowerODT
                for j=1:SampNumTemp
                    Veff(i,j)=VeffDataTemp((i-1)*SampNumTemp+j);
                end
            end
   PowerODT = 1.47:1.04125:9.8;
   TempODT = (32.6/200:32.6/200:32.6)*1e-6;
     V120091221=interp2(TempODT, PowerODT, Veff, temp, Podt,'spline');
     return;
 elseif ((Podt > 1.078) & (Podt <= 1.47)) %Interperation from lookup table 1 where the power valid range is 1.078W to 1.47W
   SampNumPowerODT=9; % There are 9 sample ODT power values: 1.078, ..., 1.47 W
   SampNumTemp=100; % There are 100 sample equlibrium temperature values: 1*1.64uK/100, 2*1.64uK/100, ... , 100*1.64uK/100
   sumdatafile=strcat(directory,'Table2V1effDec212009.dat');
   [PowerODTTemp TempODTDimpleTemp VeffDataTemp]=textread(sumdatafile, '%f%f%f','commentstyle','matlab');
   Veff=zeros(SampNumPowerODT,SampNumTemp);
            for i=1:SampNumPowerODT
                for j=1:SampNumTemp
                    Veff(i,j)=VeffDataTemp((i-1)*SampNumTemp+j);
                end
            end
   PowerODT = 1.078:0.049:1.47;
   TempODT = (1.64/100:1.64/100:1.64)*1e-6;
     V120091221=interp2(TempODT, PowerODT, Veff, temp, Podt,'spline');
     return;
else
    V120091221 =10;
    return;
end


end
