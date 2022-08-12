clear all
% The Solar-SPAM model reconstructs the solar irradiance spectrum  F, W/(m^2*nm)  
% in the 0-190 nm wavelength range at the top of the Earth's atmosphere  
% depending on the solar radio emission flux F10.7.
%
% The Solar-SPAM model was developed by Vera Nikolaeva and Evgeny Gordeev
% and published in Space Weather [Nikolaeva and Gordeev, 2022].




%% ----- INPUT -----
F107 = 123 ;  % F10.7 index within valid limits of applicability: 65 < F10.7 < 200 s.f.u. 
% -------------------------------------------------------------------------



%% ----- Solar-SPAM -----
tab_SSPAM = load('tab_SolarSPAM.txt') ;
lambda = tab_SSPAM(:,1) ;  % wavelength
p1_SSPAM = tab_SSPAM(:,2) ;
p2_SSPAM = tab_SSPAM(:,3) ;
p3_SSPAM = tab_SSPAM(:,4) ;

% -- photon energy flux at the top of atmosphere --
Ftop = p1_SSPAM.*F107^2 + p2_SSPAM.*F107 + p3_SSPAM ; 
%--------------------------------------------------------------------------


%% ----- spectrum visualization ------
figure
stem(lambda, Ftop, 'om', 'markersize', 4) ;  hold on ; grid on ;
xlim([0 190]) ; ylim([1e-6 1e-2]) ;
set(gca,'yscal','log')
xlabel('\lambda, nm') ;
ylabel('F, W m^-^2nm^-^1') ;
title('Solar-SPAM spectrum') ;
text(10, 2e-3, ['F_1_0_._7 = ',num2str(F107),' s.f.u.']) ;
%--------------------------------------------------------------------------
 



