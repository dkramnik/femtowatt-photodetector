clc
close all
clear

lw = 1.25;
fs = 14;

s = tf( 's' );

%% Noise estimate

i_ni_target = 0.42e-15;     % 0.42 fA/sqrt(Hz)

% Thorlabs FD11A
C_PD = 140e-12;
R_PD = 5e9;     % 10mV/2pA, from datasheet

C_IN = C_PD;    % Ignoring opamp Cin for now because PD is so big

R_F_vec = logspace( 0, 12, 1000 );    % Picking arb value for TIA gain for now

GBWP = 4.2e6;     % Opamp gain-BW product (took LT1792 as typical example)

C_F_vec = sqrt( C_IN ./ ( 2 * pi * GBWP * R_F_vec ) );

%disp( [ 'C_F = ' num2str( C_F * 1e12 ) ' pF' ] );

f_3dB_TIA_vec = sqrt( GBWP ./ ( 2 * pi * R_F_vec * C_IN ) );

%disp( [ 'f_3dB_TIA = ' num2str( f_3dB_TIA * 1e-3 ) ' kHz' ] );

% I have adjusted the opamp noise terms so that they never dominate,
% regardless of R_F selection
% These are the specs we should try to meet in the opamp selection

i_n_opamp = 1.5e-15;    % 1.5 fA/sqrt(Hz)
e_n_opamp = 1.0e-9;     % 1.0 nV/sqrt(Hz)
k_B = physconst( 'Boltzmann' );
T = 300;

i_ni_pd_resistance_term = sqrt( 4 * k_B * T ./ R_PD ) * ones( size( R_F_vec ) );
i_ni_opamp_i_term = i_n_opamp * ones( size( R_F_vec ) );
i_ni_opamp_e_term = ( e_n_opamp ./ R_F_vec );
i_ni_R_F_term = sqrt( 4 * k_B * T ./ R_F_vec );
i_ni_peaking_term = sqrt( ( e_n_opamp * 2 * pi * f_3dB_TIA_vec * C_IN ).^2 / 3 );

i_ni_vec = sqrt( ...
    i_ni_pd_resistance_term.^2 + ...
    i_ni_opamp_i_term.^2 + ...
    i_ni_opamp_e_term.^2 + ...
    i_ni_R_F_term.^2 + ...
    i_ni_peaking_term.^2 ...
    );

%disp( [ 'Input-referred noise = ' num2str( i_ni * 1e15 ) ' fA/sqrt(Hz)' ] );

fig1 = figure( );

subplot( 2, 1, 1 ); hold all; grid on;

plot( R_F_vec, 1e15 * i_ni_vec, 'linewidth', 2 );

plot( R_F_vec, 1e15 * i_ni_pd_resistance_term, '--', 'linewidth', lw );
plot( R_F_vec, 1e15 * i_ni_R_F_term, '--', 'linewidth', lw );
plot( R_F_vec, 1e15 * i_ni_opamp_i_term, '--', 'linewidth', lw );
plot( R_F_vec, 1e15 * i_ni_opamp_e_term, '--', 'linewidth', lw );
plot( R_F_vec, 1e15 * i_ni_peaking_term, '--', 'linewidth', lw );

plot( R_F_vec, 1e15 * i_ni_target * ones( size( R_F_vec ) ), 'k', 'linewidth', lw );

legend( ...
    'Total', ...
    'R_{PD} Term (Thorlabs FD11A)', ...
    'R_F Term', ...
    'Opamp Current Term', ...
    'Opamp Voltage Term', ...
    'C_{IN} Peaking Term', ...
    'Target' ...
    );

set( gca, 'xscale', 'log' );
set( gca, 'yscale', 'log' );

ylim( [ 1e-1, 1e6 ] );

xlabel( 'R_F' );
ylabel( 'Input-Referred Noise [fA/sqrt(Hz)]' );

set( gca, 'fontsize', fs );

subplot( 2, 1, 2 ); hold all; grid on;

plot( R_F_vec, f_3dB_TIA_vec, 'linewidth', lw );

set( gca, 'xscale', 'log' );
set( gca, 'yscale', 'log' );

xlabel( 'R_F' );
ylabel( 'Ideal TIA BW [Hz]' );

set( gca, 'fontsize', fs );

%% Loop Shaping

fp1 = 10;
Aol = 1.3e6 / ( s / ( 2 * pi * fp1 ) + 1 );   % 100MHz GBP

Z_i = prl( R_PD, 1 / ( s * C_PD ) );
Z_f = 100e9;
fb_div = Z_i / ( Z_i + Z_f );

Tloop = fb_div * Aol;
Acl = prl( Z_i, Z_f ) * Aol / ( 1 + Aol * Z_i / ( Z_i + Z_f ) );

fig2 = figure( ); hold all;

opts = bodeoptions;
opts.FreqUnits = 'Hz';

bode( Aol, opts );
bode( Tloop, opts );
bode( Acl, opts );

legend( 'Opamp Open-Loop Gain', 'TIA Loop Gain', 'TIA Closed-Loop Transimpedance' );
grid on;


