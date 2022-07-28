function [qN,qN2,qO,qO2,qNO] = fun_StraightPhotoionization(INPUT)
% fun_StraightPhotoionization.m ::: calculates the photoionization rates Q of 
% the basic atmospheric gas molecules and atoms (O, O2, N, N2, NO) 
% due to impact of the straigt solar radiation.
%
% More details on the photoionization calculation procedure can be found in
% [Nikolaeva et al., Atmosphere 2021], and on the SPAM solar spectrum model in
% [Nikolaeva & Gordeev, GMD 2022]
%
% INPUT is a structure array containing all needed input variables


%% ----- INPUT -----
alts = INPUT.alts ;  % altitude range for Q calculation

% -- MSIS'00 in altitude range [79:1000]km --
Year = INPUT.Year ; 
Mon = INPUT.Mon ;
Day = INPUT.Day ;
Hour = INPUT.Hour ;  
alt_msis = INPUT.alt_msis ;
GLat = INPUT.GLat ;  % Geodetic Latitude range -90:90
GLon = INPUT.GLon ;  % Geodetic Longitude range -180:+180
nO_msis = INPUT.nO_msis ;   % m-3
nN2_msis = INPUT.nN2_msis ;
nO2_msis = INPUT.nO2_msis ;
T_msis = INPUT.T_msis ;
nHe_msis = INPUT.nHe_msis ; 
nAr_msis = INPUT.nAr_msis ;
nH_msis = INPUT.nH_msis ;
nN_msis = INPUT.nN_msis ;
F107 = INPUT.F107 ;

% -- NO empirical concentration by [Mitra, 1968] --
nNO_mitra = INPUT.nNO_mitra ;  
% -------------------------------------------------------------------------



%% ----- some universal constants -----
kB = 1.3805404e-23 ;  % Boltzmann const
amu = 1.660539e-27 ;  % Unified atomic mass unit
%--------------------------------------------------------------------------



%% ----- Aero-SPAM [Nikolaeva and Gordeev, 2022] -----
tab_ASPAM = load('tab_AeroSPAM.txt') ;
p1_ASPAM = tab_ASPAM(:,4) ;
p2_ASPAM = tab_ASPAM(:,5) ;
p3_ASPAM = tab_ASPAM(:,6) ;

% -- photon flux at the top of atmosphere --
Itop_UV = p1_ASPAM(1:end-1).*F107^2 + p2_ASPAM(1:end-1).*F107 + p3_ASPAM(1:end-1) ; 
Itop_La = p1_ASPAM(end).*F107^2 + p2_ASPAM(end).*F107 + p3_ASPAM(end) ;  % Ly-alpha flux
%--------------------------------------------------------------------------



%% ----- Ionization Cross Sections [Nikolaeva et al., 2021, Appendix C] -----
tab_CS = load('tab_CrossSections.txt') ;
sO_i = tab_CS(:,3) *1e-22 ;  % sigma O -- cross-sect  i -- ionization    [m2]  
sO_a = tab_CS(:,4) *1e-22 ;   % sigma O -- cross-sect  a -- absorbtion   [m2]
sN_i = tab_CS(:,5) *1e-22 ;
sN_a = tab_CS(:,6) *1e-22 ;
sO2_i = tab_CS(:,7) *1e-22 ;
sO2_a = tab_CS(:,8) *1e-22 ;
sN2_i = tab_CS(:,9) *1e-22 ;
sN2_a = tab_CS(:,10) *1e-22 ;
sNO_i = tab_CS(:,11) *1e-22 ;
sNO_a = tab_CS(:,12) *1e-22 ;
%--------------------------------------------------------------------------



%% ----- Chapman function calculation -----
DoY = datenum(Year,Mon,Day)-datenum(Year,0,0) ; % Day of Year
declin = -23.44*cosd((360/365)*(DoY+10)) ; % solar declination   
lst = Hour + GLon/15 ;  % Local apparent solar time [nrmlsis-2000]
sha = 15*(lst-12) ; % solar hour angle
SZA = 90-asind( sind(GLat)*sind(declin) + cosd(GLat)*cosd(declin)*cosd(sha)) ; % Solar Zenith Angle

Chap = zeros(size(alt_msis)) ;  % Chapman Function
for k = 1:length(alt_msis) 
    MolWeights = [4.0026, 15.9994, 28.0134, 31.9988, 39.948, 1.00797, 14.0067, 30.0061] ;  % Molecular Weights [g/mol]
    NeutralsDens = [nHe_msis(k), nO_msis(k), nN2_msis(k), nO2_msis(k), nAr_msis(k), nH_msis(k), nN_msis(k), nNO_mitra(k)] ;  % concentration of neutrals 
    sma = sum(NeutralsDens.*MolWeights)/sum(NeutralsDens) ; % specific mass of the air [g/mol]
    grav = 9.780327*(1+0.0053024*sind(GLat)^2-0.0000058*sind(2*GLat)^2)-3.086*1e-6*alt_msis(k) ; % gravitational acceleration
    ScH = (kB*T_msis(k))/(sma*amu*grav) ;  % scale height
    XX = (6371+alt_msis(k))*1e3/ScH ;  
    AEF = sqrt(XX/2) * abs(cosd(SZA)) ; % Additional Error Function
    EF = (1.0606963+0.55643831*AEF) / (1.0619896 + 1.7245609*AEF + AEF^2) * (exp(1)^(-AEF^2)) ;  % Error Function
    if SZA <= 90 
        Chap(k) = ((pi/2*XX)^0.5) * (exp(1)^(AEF^2)) * EF ;  % Chapman Function
    elseif SZA > 90
        Chap(k) = ((2*pi*XX)^0.5) * (sind(SZA)^0.5*exp(XX*(1-sind(SZA)))-0.5*exp(1)^(AEF^2)*EF) ;  
    end  
end
%--------------------------------------------------------------------------



%% ----- photoionization rates vertical profile -----
for ka = 1:length(alts) 
    
    % -- values at a specified height --
    f = find(alt_msis==alts(ka)) ; % 
    alt = alt_msis(f) ;
    nO = nO_msis(f) ;   % [m-3]
    nN2 = nN2_msis(f) ;
    nO2 = nO2_msis(f) ;
    nN = nN_msis(f) ;
    nNO = nNO_mitra(f) ;
    Talt = T_msis(f) ;   % Atm temp at altititude

    % -- Height integrated concentration --
    f = find(alt_msis>=alt) ;
    integr_nO = sum(Chap(f).*nO_msis(f).*(alt_msis(f)-alt_msis(f-1)).*1e3) ;  % integral concentration in the above atmospheric column
    integr_nO2 = sum(Chap(f).*nO2_msis(f).*(alt_msis(f)-alt_msis(f-1)).*1e3) ;
    integr_nN2 = sum(Chap(f).*nN2_msis(f).*(alt_msis(f)-alt_msis(f-1)).*1e3) ;
    integr_nN = sum(Chap(f).*nN_msis(f).*(alt_msis(f)-alt_msis(f-1)).*1e3) ;
    integr_nNO = sum(Chap(f).*nNO_mitra(f).*(alt_msis(f)-alt_msis(f-1)).*1e3) ;
    clear f ;

    % -- photoionization rates, m-3s-1 --
    qO(ka) = nO * sum(sO_i .* Itop_UV .* exp(-(sO_a*integr_nO + sO2_a*integr_nO2 + sN2_a*integr_nN2 + sN_a*integr_nN + sNO_a*integr_nNO) )) ;
    qO2(ka) = nO2 * sum(sO2_i .* Itop_UV .* exp(-(sO_a*integr_nO + sO2_a*integr_nO2 + sN2_a*integr_nN2 + sN_a*integr_nN + sNO_a*integr_nNO) )) ;
    qN(ka) = nN * sum(sN_i .* Itop_UV .* exp(-(sO_a*integr_nO + sO2_a*integr_nO2 + sN2_a*integr_nN2 + sN_a*integr_nN + sNO_a*integr_nNO) )) ;
    qN2(ka) = nN2 * sum(sN2_i .* Itop_UV .* exp(-(sO_a*integr_nO + sO2_a*integr_nO2 + sN2_a*integr_nN2 + sN_a*integr_nN + sNO_a*integr_nNO) )) ;
    qNO_UV(ka) = nNO * sum(sNO_i .* Itop_UV .* exp(-(sO_a*integr_nO + sO2_a*integr_nO2 + sN2_a*integr_nN2 + sN_a*integr_nN + sNO_a*integr_nNO) )) ;
    qNO_La(ka) = nNO * sum(1.944e-22 .* Itop_La .* exp(-(9e-25*integr_nO2 + 6e-27*integr_nN2 + 2.4e-22*integr_nNO) )) ;
    qNO(ka) = qNO_UV(ka) + qNO_La(ka) ;
end
%--------------------------------------------------------------------------

