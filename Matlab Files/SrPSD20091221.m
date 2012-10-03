function PSDCal20091221=SrPSD20091221(Podt,temp,NumAtoms,MassAtom)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSDCal20091221 = SrPSD20091221(PowerODT, temp, NumAtoms, MassAtom) for
%% Dec. 21, 2009 ODT setup
%% where Powers are in Watts, temp in Kelvin and Mass in Atomic Weight e.g. input 88 for 88Sr.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h=6.626*10^(-34); %Planck's constant
kBoltz=1.38*10^(-23); %Blotzmann constant
mass=MassAtom.*1.672.*10.^(-27); %atom mass
lambdadB=h./sqrt(2.*pi.*mass.*kBoltz.*temp); %de Broglie wavelength
PSDCal20091221=NumAtoms./V1eff20091221(Podt,temp).*lambdadB.^3;
        
end

