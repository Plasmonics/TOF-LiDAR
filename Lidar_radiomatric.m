%%    Lidar_Radiometric_analysis

%% ------------------------------------
%  TOF LiDAR system signal processing for static run
%  Last update: 11/16/2020
%  Author:Hosna Sultana
% For TOF LiDAR timestamp data (time and number of count) for getting range from the maximum probability
% density from the histogram and multiple target range from most prominent peaks

%% ------------------------------------

Pp = 100;   % power of the laser doide in a pulse in W
tau = 10 * 10^-9 ; f_rep = 1 *10^3 ;              % pulse width  and rep frequency
P_ave = Pp * tau * f_rep ;                        % average power
T_L_asp = 0.95; T_20x_BE = 0.94; T_3x_BE = 0.94;  T_filt = 0.90;    % all transmission coefficients
Po_p = Pp * T_L_asp * T_20x_BE;                   % power of the outgoing signal in a pulse in W
Po_ave = P_ave * T_L_asp * T_20x_BE;              % Average power of the outgoing signal in a pulse in W
M = 1;                                % atmospheric absorption coefficient
rho = 0.55;                         % reflectivity of the target(frozen snow with grain size 100 micron at 905 nm wavelength)
AT = 0.046 * 0.046;                 % area of target in m^2
At = 0.046 * 0.046;                 % Foot print illuminated area of target in m^2
Ad = pi * ((180*10^-6)/2)^2 ;       % area of the detector in m^2
Ar = pi * ((23*10^-3)/2)^2;         % area of the receiver in m^2
Af = pi * ((8*10^-3)/2)^2;          % area of the receiver in m^2
R = 500 ;                           % range 500 m
 
Phi_g = Po_p* M ;                     % power on fall on ground (W)
phi_t = Po_p * M / At  ;              %(W/m2) Irradiance at target 
phi_pt = phi_t * AT  ;              %Power of the transmitted pulse on the target:  ?t AT = [Po M/ At] AT (W)
phi_rf_t = phi_pt * rho;            %Power of the reflected pulse from the target:  ?t AT ? = [Po M/ At] AT  ?(W)
% phi_rf_t_s = phi_rf_t *M / pi
phi_rf_t_s = phi_rf_t *M / pi;      %Power of the reflected pulse received at receiver from per unit solid angle the 
                                    %target:?t AT ? M/ ? = [Po M2/ At] AT ? /?  (W)
phi_r_r_s = (Ar/ R^2)*  phi_rf_t_s ;  %  Power received at receiver in complete solid angle of receiver:  
                                      % [Ar/ R2 ]?t AT ? M2/ At ? = Po M2AT ? Ar / At ?  R2      (W/m2)
                                                                                                                 
phi_r_D =  phi_r_r_s *  T_3x_BE *  T_filt;   %Power received at detector in complete solid angle of receiver (W)
phi_r_D_active = phi_r_D * (Ad/Ar) ;       %fraction of Power received at detector active area in complete solid angle fromreceiver (W)

%% DETECTOR

          % Number of photon/ sec = 5.03 * 10^15 * lanbda in nm * optical power (W)
          % Number of photon/32 nano sec = [5.03 * 10^15 * lanbda in nm * optical power (W)]* (32 / 10^9)
          
N_p = 5.03 * 10^15 * 905 * phi_r_D_active          
%N_p = 5.03 * 10^15 * 905 * phi_r_D_active * (32 / 10^9)

t_d = 22* 10^-9;               % (s) Module dead time
c_r = 37 *10^6 ;               %(c/s)output count rate
d_C_r = 500 ;                  % (c/s) or (cps) dark count rate 
PDE = 0.40 ;                   % at 905 nm, photon detection efficiency, from graph 4 of the photon detector manual
C_F = 1/( 1-(t_d * c_r))  ;    % Correction factor
M_C_r = 8 *10^6 ;            % convert N_p to Module count rate from graph 7 of the photon detector manual

%% actual count rate = 
       % [(output mudule count rate * correction factor @ module count rate) - Dark count rate]/ photon detection efficiency module

 a_C_r = ((M_C_r * C_F)- d_C_r) / PDE      % Acctual photon count rate

 










