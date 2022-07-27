clear all
% The Aero-SPAM model reconstructs the shortwave photon flux I, 1/(m^2*s*nm)  
% at the top of the Earth's atmosphere depending on the solar radio emission 
% flux F10.7.
%
% The Aero-SPAM model was developed by Vera Nikolaeva and Evgeny Gordeev
% and published in Geoscientific Model Development (GMD)
% [Nikolaeva and Gordeev, 2022].
%
% Note that the visualization is performed for specific spectral channels 
% numbered from 1 to 37, which include photon flux values in 17 spectral 
% lines and 20 spectral bands (see tab_AeroSPAM.txt).



%% ----- INPUT -----
F107 = 200 ;  % F10.7 index within valid limits of applicability: 65 < F10.7 < 200 s.f.u. 
% -------------------------------------------------------------------------



%% ----- Aero-SPAM -----
tab_ASPAM = load('tab_AeroSPAM.txt') ;
chan = tab_ASPAM(:,1) ;  % spectral channel
p1_ASPAM = tab_ASPAM(:,4) ;
p2_ASPAM = tab_ASPAM(:,5) ;
p3_ASPAM = tab_ASPAM(:,6) ;

% -- photon flux at the top of atmosphere --
Itop = p1_ASPAM.*F107^2 + p2_ASPAM.*F107 + p3_ASPAM ; 
%--------------------------------------------------------------------------


%% ----- spectrum visualization ------
figure
stem(chan, Itop, 'or') ;  hold on ; grid on ;
xlim([0 38]) ; ylim([1e10 1e16]) ;
set(gca,'yscal','log')
xlabel('Channel #') ;
ylabel('I, m^-^2s^-^1nm^-^1') ;
title('Aero-SPAM "spectrum"') ;
text(2, 2e+15, ['F_1_0_._7 = ',num2str(F107),' s.f.u.']) ;
%--------------------------------------------------------------------------
 



