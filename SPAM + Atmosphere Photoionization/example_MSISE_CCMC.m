clear all

% Calculation of photoionization rate profiles with the NRLMSISE-00
% atmosphere results obtained from the Comunity Coordinated Modelling Center 
% webpage https://ccmc.gsfc.nasa.gov/

% Format of the NRLMSISE-00 file:
% 1) Altitude range 79:1:1000 km;
% 2) Required fields (select when prepare the MSISE output at CCMC) are specified 
% below in the NPUT section. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ----- INPUT -----
alts = 80:250 ;  % altitude range for Q calculation

% -- MSIS'00 in altitude range [79:1000]km --
% msis = load('nrlmsise_2009Jun20_10UT.txt') ;  % neutral atmosphere model(F10.7=65 sfu)
msis = load('nrlmsise_2015Jun18_10UT.txt') ;  % neutral atmosphere model(F10.7=155 sfu)
Year = msis(1,1) ;               % Year
Mon = msis(1,2) ;                % Month
Day = msis(1,3) ;                % Day (of Month)
Hour = msis(1,4) ;               % Hour (of Day, decimal)
alt_msis = msis(:,5) ;           % Height, km
GLat = msis(1,6) ;               % Geographic Latitude range -90:90
GLon = msis(1,7) ;               % Geographic Longitude range -180:+180
nO_msis = msis(:,8) .*1e+6 ;     % Concentration of Atomic Oxygen O [cm-3 --> m-3]
nN2_msis = msis(:,9) .*1e+6 ;    % Nitrogen N2
nO2_msis = msis(:,10) .*1e+6 ;   % Oxygen O2
T_msis = msis(:,11) ;            % Neutral Tempereature [K]
nHe_msis = msis(:,12) .*1e+6 ;   % Helium He
nAr_msis = msis(:,13) .*1e+6 ;   % Argon Ar
nH_msis = msis(:,14) .*1e+6 ;    % Hydrogen H
nN_msis = msis(:,15) .*1e+6 ;    % Nitrogen N
F107 = msis(1,16) ;              % F10.7 (daily)

% -- NO empirical concentration by [Mirta, ...] --
nNO_mirta = 0.4*exp(-3700./T_msis).*nO2_msis + 5*1e-7*nO_msis ;  

% -- combined Input --
INPUT.alts = alts ;      INPUT.Year = Year ;        INPUT.Mon = Mon ;        INPUT.Day = Day ;
INPUT.Hour = Hour ;      INPUT.alt_msis = alt_msis ;    INPUT.GLat = GLat ;       INPUT.GLon = GLon ;
INPUT.nO_msis = nO_msis ;   INPUT.nN2_msis = nN2_msis ;   INPUT.nO2_msis = nO2_msis ;  INPUT.T_msis = T_msis ;
INPUT.nHe_msis = nHe_msis ; INPUT.nAr_msis = nAr_msis ;   INPUT.nH_msis = nH_msis ;   INPUT.nN_msis = nN_msis ;
INPUT.F107 = F107 ;     INPUT.nNO_mirta = nNO_mirta ;
% -------------------------------------------------------------------------


%% ----- call photoionization function -----
[qN,qN2,qO,qO2,qNO] = fun_StraightPhotoionization(INPUT) ;
% -------------------------------------------------------------------------


%% ----- plot results in [cm-3 s-1] -----
figure
semilogx(qN.*1e-6,alts,'c') ; hold on ; grid on ; % qN
semilogx(qN2.*1e-6,alts,'r') ;  % qN2
semilogx(qO.*1e-6,alts,'b') ;  % qO
semilogx(qO2.*1e-6,alts,'m') ;  % qO2
semilogx(qNO.*1e-6,alts,'g') ;  % qNO
semilogx(sum([qN;qN2;qO;qO2;qNO]).*1e-6,alts,'k', 'linewidth', 2) ;  % q total
xlim([1e+0 1e+4]) ; 
ylim([alts(1) alts(end)]) ;
text(15,max(alts)-15,'qN', 'color','c') ;
text(15,max(alts)-25,'qN_2', 'color','r') ;
text(15,max(alts)-35,'qO', 'color','b') ;
text(15,max(alts)-45,'qO_2', 'color','m') ;
text(15,max(alts)-55,'qNO', 'color','g') ;
text(15,max(alts)-65,'\Sigmaq', 'color','k') ;
xlabel('q, cm^-^3s^-^1') ;
ylabel('Altitude, km') ;

