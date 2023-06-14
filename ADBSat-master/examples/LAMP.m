% Cube Example
%
%--- Copyright notice ---%
% Copyright (C) 2021 The University of Manchester
% Written by David Mostaza Prieto,  Nicholas H. Crisp, Luciana Sinpetru and Sabrina Livadiotti
%
% This file is part of the ADBSat toolkit.
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%------------ BEGIN CODE ----------%

clear


%Input conditions
alt = 300; %km

% Model parameters
inparam.gsi_model = 'sentman';
%inparam.alpha = 1; % Accommodation (altitude dependent)
inparam.Tw = 300; % Wall Temperature [K] -> jose

inparam.alpha = 0.9; % Accommodation (altitude dependent), 0.9 is more limiting from calculations
inc = 96.67; %deg
% 238 367 average for 65 65 

% for lat = -90:1:90 
%     for lon = -90:1:90
env = [alt*1e3, inc/2, 0, 106, 0, 238, 367, ones(1,7)*3, 0]; % Environment variables
% Environment Calculations
inparam = environment(inparam, env(1),env(2),env(3),env(4),env(5),env(6),env(7),env(8:14),env(15));

cg_s = {[0, 0.181, -0.018], [0, 0.177, -0.018], [0, 0.172, -0.018]};  % array of cg's y, z, x for blender its weird how it works

%cg_s = {[0, 0.177, 0.]};  % array of cg's y, z, x for blender its weird how it works

% cg_s = {};
% for i = 0:0.002:0.01
%     for j = 0:0.002:0.01:
%         for k = 0:0.002:0.01:
%             for cg = cg_nominal:
% 
%             cg_s = [cg_s; [i, j, k]];
%         end
%     end
% end


modGeneral = 'CI-Br-Sf0';
modNames = {'CI-Br0-Sf0', 'CI-Br30-Sf0', 'CI-Br60-Sf0', 'CI-Br90-Sf0'};
config_1_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CI-Br-Sf0\\CI-Br-Sf0.csv';
get_config_torque_csv(config_1_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CIII-Bfhi-Sfi';
modNames = {'CIII-Bfh0-Sf0', 'CIII-Bfh12-Sf12'};
config_2_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CIII-Bfhi-Sfi\\CIII-Bfhi-Sfi.csv';
get_config_torque_csv(config_2_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CIII-Bfhi-Sfo';
modNames = {'CIII-Bfh0-Sfo', 'CIII-Bfh12-Sfo'};
config_3_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CIII-Bfhi-Sfo\\CIII-Bfhi-Sfo.csv';
get_config_torque_csv(config_3_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CIII-Bfvi-Sfi';
modNames = {'CIII-Bfv46-Sf46', 'CIII-Bfv61-Sf61', 'CIII-Bfv72-Sf72', 'CIII-Bfv80-Sf80', 'CIII-Bfv83-Sf83'};
config_4_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CIII-Bfvi-Sfi\\CIII-Bfvi-Sfi.csv';
get_config_torque_csv(config_4_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CIII-Bfvi-Sf0';
modNames = {'CIII-Bfv46-Sf0', 'CIII-Bfv61-Sf0', 'CIII-Bfv72-Sf0', 'CIII-Bfv80-Sf0', 'CIII-Bfv83-Sf0'};
config_5_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CIII-Bfvi-Sf0\\CIII-Bfvi-Sf0.csv';
get_config_torque_csv(config_5_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CIV-Ffhi-Sfi';
modNames = {'CIV-Ffh0-Sf0', 'CIV-Ffh12-Sf12'};
config_6_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CIV-Ffhi-Sfi\\CIV-Ffhi-Sfi.csv';
get_config_torque_csv(config_6_loc, cg_s, inparam, modNames, modGeneral)

% CI-Bri-Sfi
% CI-Br-Sf12
modGeneral = 'CI-Br-Sfi\\CI-Br-Sf12';
modNames = {'CI-Br0-Sf12', 'CI-Br30-Sf12', 'CI-Br60-Sf12', 'CI-Br90-Sf12'};
config_7_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CI-Br-Sfi\\CI-Br-Sf12\\CIV-Br-Sf12.csv';
get_config_torque_csv(config_7_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CI-Br-Sfi\\CI-Br-Sf46';
modNames = {'CI-Br0-Sf46', 'CI-Br30-Sf46', 'CI-Br60-Sf46', 'CI-Br90-Sf46'};
config_8_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CI-Br-Sfi\\CI-Br-Sf46\\CIV-Br-Sf46.csv';
get_config_torque_csv(config_8_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CI-Br-Sfi\\CI-Br-Sf61';
modNames = {'CI-Br0-Sf61', 'CI-Br30-Sf61', 'CI-Br60-Sf61', 'CI-Br90-Sf61'};
config_9_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CI-Br-Sfi\\CI-Br-Sf61\\CIV-Br-Sf61.csv';
get_config_torque_csv(config_9_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CI-Br-Sfi\\CI-Br-Sf71';
modNames = {'CI-Br0-Sf72', 'CI-Br30-Sf72', 'CI-Br60-Sf72', 'CI-Br90-Sf72'};
config_10_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CI-Br-Sfi\\CI-Br-Sf72\\CIV-Br-Sf72.csv';
get_config_torque_csv(config_10_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CI-Br-Sfi\\CI-Br-Sf80';
modNames = {'CI-Br0-Sf80', 'CI-Br30-Sf80', 'CI-Br60-Sf80', 'CI-Br90-Sf80'};
config_11_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CI-Br-Sfi\\CI-Br-Sf80\\CIV-Br-Sf80.csv';
get_config_torque_csv(config_11_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CI-Br-Sfi\\CI-Br-Sf83';
modNames = {'CI-Br0-Sf83', 'CI-Br30-Sf83', 'CI-Br60-Sf83', 'CI-Br90-Sf83'};
config_12_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CI-Br-Sfi\\CI-Br-Sf83\\CIV-Br-Sf83.csv';
get_config_torque_csv(config_12_loc, cg_s, inparam, modNames, modGeneral)


% CIV-Fr-Sfi
modGeneral = 'CIV-Fr-Sfi\\CIV-Fr-Sf0';
modNames = {'CIV-Fr0-Sf0', 'CIV-Fr30-Sf0', 'CIV-Fr60-Sf0', 'CIV-Fr90-Sf0'};
config_13_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CIV-Fr-Sfi\\CIV-Fr-Sf0\\CIV-Ffhi-Sf0.csv';
get_config_torque_csv(config_13_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CIV-Fr-Sfi\\CIV-Fr-Sf12';
modNames = {'CIV-Fr-Sf12', 'CIV-Fr30-Sf12', 'CIV-Fr60-Sf12', 'CIV-Fr90-Sf12'};
config_14_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CIV-Fr-Sfi\\CIV-Fr-Sf12\\CIV-Ffhi-Sf12.csv';
get_config_torque_csv(config_14_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CIV-Fr-Sfi\\CIV-Fr-Sf46';
modNames = {'CIV-Fr0-Sf46', 'CIV-Fr30-Sf46', 'CIV-Fr60-Sf46', 'CIV-Fr90-Sf46'};
config_15_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CIV-Fr-Sfi\\CIV-Fr-Sf46\\CIV-Ffhi-Sf46.csv';
get_config_torque_csv(config_15_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CIV-Fr-Sfi\\CIV-Fr-Sf61';
modNames = {'CIV-Fr0-Sf61', 'CIV-Fr30-Sf61', 'CIV-Fr60-Sf61', 'CIV-Fr90-Sf61'};
config_16_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CIV-Fr-Sfi\\CIV-Fr-Sf61\\CIV-Fr-Sf61.csv';
get_config_torque_csv(config_16_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CIV-Fr-Sfi\\CIV-Fr-Sf71';
modNames = {'CIV-Fr0-Sf72', 'CIV-Fr30-Sf72', 'CIV-Fr60-Sf72', 'CIV-Fr90-Sf72'};
config_17_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CIV-Fr-Sfi\\CIV-Fr-Sf72\\CIV-Fr-Sf72.csv';
get_config_torque_csv(config_17_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CIV-Fr-Sfi\\CIV-Fr-Sf80';
modNames = {'CIV-Fr0-Sf80', 'CIV-Fr30-Sf80', 'CIV-Fr60-Sf80', 'CIV-Fr90-Sf80'};
config_18_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CIV-Fr-Sfi\\CIV-Fr-Sf80\\CIV-Fr-Sf80.csv';
get_config_torque_csv(config_18_loc, cg_s, inparam, modNames, modGeneral)

modGeneral = 'CIV-Fr-Sfi\\CIV-Fr-Sf83';
modNames = {'CIV-Fr0-Sf83', 'CIV-Fr30-Sf83', 'CIV-Fr60-Sf83', 'CIV-Fr90-Sf83'};
config_19_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\obj_files\\Configs\\CIV-Fr-Sfi\\CIV-Fr-Sf83\\CIV-Fr-Sf83.csv';
get_config_torque_csv(config_19_loc, cg_s, inparam, modNames, modGeneral)



function get_config_torque_csv(config_loc, cg_s, inparam, modNames, modGeneral)
    ADBSat_path = ADBSat_dynpath;
    config_f_id = fopen(config_loc,'w');
    fprintf(config_f_id, 'case, max torque roll, max torque pitch, max torque yaw, max unstable torque pitch, max unstable torque yaw, max torque roll prop, max torque pitch prop, max torque yaw prop, max unstable torque pitch prop, max unstable torque yaw prop, max average roll torque all, max average pitch torque all, max average yaw torque all, max roll torque all, max pitch torque all, max yaw torque all\n');
    [abs_t_m_r, abs_t_m_p, abs_t_m_y, abs_t_m_p_u, abs_t_m_y_u, abs_t_m_r_prop, abs_t_m_p_prop, abs_t_m_y_prop, abs_t_m_p_u_prop, abs_t_m_y_u_prop, abs_t_m_avg_r, abs_t_m_avg_p, abs_t_m_avg_y, abs_t_m_all_r, abs_t_m_all_p, abs_t_m_all_y] = configs_torque_calculator(ADBSat_path, modNames, modGeneral, cg_s, inparam, config_f_id);
    fprintf(config_f_id, 'absolute maximums\n');
    formatEnd = '%8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f \n';
    fprintf(config_f_id, formatEnd, abs_t_m_r, abs_t_m_p, abs_t_m_y, abs_t_m_p_u, abs_t_m_y_u, abs_t_m_r_prop, abs_t_m_p_prop, abs_t_m_y_prop, abs_t_m_p_u_prop, abs_t_m_y_u_prop, abs_t_m_avg_r, abs_t_m_avg_p, abs_t_m_avg_y, abs_t_m_all_r, abs_t_m_all_p, abs_t_m_all_y);
    fclose(config_f_id);
end

function [abs_t_m_r, abs_t_m_p, abs_t_m_y, abs_t_m_p_u, abs_t_m_y_u, abs_t_m_r_prop, abs_t_m_p_prop, abs_t_m_y_prop, abs_t_m_p_u_prop, abs_t_m_y_u_prop, abs_t_m_avg_r, abs_t_m_avg_p, abs_t_m_avg_y, abs_t_m_all_r, abs_t_m_all_p, abs_t_m_all_y] = configs_torque_calculator(ADBSat_path, modNames, modGeneral, cg_s, inparam, config_1_f_id)
    abs_t_m_r = 0; % absolute torque max roll
    abs_t_m_p = 0; % absolute torque max pitch
    abs_t_m_y = 0; % absolute torque max yaw
    
    abs_t_m_p_u = 0; % absolute torque max pitch unstable
    abs_t_m_y_u = 0; % absolute torque max yaw unstable

    abs_t_m_r_prop = 0; % absolute torque max roll for propulsion on
    abs_t_m_p_prop = 0; % absolute torque max pitch for propulsion on
    abs_t_m_y_prop = 0; % absolute torque max yaw for propulsion on
    
    abs_t_m_p_u_prop = 0; % absolute torque max pitch unstable for propulsion on
    abs_t_m_y_u_prop = 0; % absolute torque max yaw unstable for propulsion on

    abs_t_m_avg_r = 0;
    abs_t_m_avg_p = 0;
    abs_t_m_avg_y = 0;
    abs_t_m_all_r = 0;
    abs_t_m_all_p = 0;
    abs_t_m_all_y = 0;

    formatRes = '%s, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f \n';

    for modName_el = modNames
        modName = char(modName_el);
        display(modName)
        modIn = fullfile(ADBSat_path,'inou','obj_files','Configs', modGeneral,[modName,'.obj']);
        modOut = fullfile(ADBSat_path,'inou','models');
        resOut = fullfile(ADBSat_path,'inou','results', modName);
    
        [t_m_r, t_m_p, t_m_y, t_m_p_u, t_m_y_u, t_m_r_prop, t_m_p_prop, t_m_y_prop, t_m_p_u_prop, t_m_y_u_prop, t_m_avg_r, t_m_avg_p, t_m_avg_y, t_m_all_r, t_m_all_p, t_m_all_y] = calc_adcs_torques(cg_s, inparam, modIn, modOut, resOut);
        % display(t_m_r)
        % display(t_m_p)
        % display(t_m_y)
        % display(t_m_p_u)
        % display(t_m_y_u)
    
        fprintf(config_1_f_id, formatRes, modName, t_m_r, t_m_p, t_m_y, t_m_p_u, t_m_y_u, t_m_r_prop, t_m_p_prop, t_m_y_prop, t_m_p_u_prop, t_m_y_u_prop,  t_m_avg_r, t_m_avg_p, t_m_avg_y, t_m_all_r, t_m_all_p, t_m_all_y);

        abs_t_m_r = max(abs_t_m_r, t_m_r);
        abs_t_m_p = max(abs_t_m_p, t_m_p); % 3!!!!!!! around z here because either wrong model or wrong program
        abs_t_m_y = max(abs_t_m_y, t_m_y); 
            
        abs_t_m_p_u = max(abs_t_m_p_u, t_m_p_u);
        abs_t_m_y_u = max(abs_t_m_y_u, t_m_y_u);

        abs_t_m_r_prop = max(abs_t_m_r_prop, t_m_r_prop);
        abs_t_m_p_prop = max(abs_t_m_p_prop, t_m_p_prop); % 3!!!!!!! around z here because either wrong model or wrong program
        abs_t_m_y_prop = max(abs_t_m_y_prop, t_m_y_prop); 
            
        abs_t_m_p_u_prop = max(abs_t_m_p_u_prop, t_m_p_u_prop);
        abs_t_m_y_u_prop = max(abs_t_m_y_u_prop, t_m_y_u_prop);

        abs_t_m_avg_r = max(abs_t_m_avg_r, t_m_avg_r);
        abs_t_m_avg_p = max(abs_t_m_avg_p, t_m_avg_p);
        abs_t_m_avg_y = max(abs_t_m_avg_y, t_m_avg_p);

        abs_t_m_all_r = max(abs_t_m_all_r, t_m_avg_r);
        abs_t_m_all_p = max(abs_t_m_all_p, t_m_avg_p);
        abs_t_m_all_y = max(abs_t_m_all_y, t_m_avg_y);

    end
end

function [max_roll, max_pitch, max_yaw, max_pitch_unstable, max_yaw_unstable, max_p_roll, max_p_pitch, max_p_yaw, max_p_pitch_unstable, max_p_yaw_unstable, max_avg_roll, max_avg_pitch, max_avg_yaw, max_all_roll, max_all_pitch, max_all_yaw] = calc_adcs_torques(cg_s, inparam, modIn, modOut, resOut)
    
    shadow = 1;

    solar = 0;
    %inparam.sol_cR = 0.77; % Specular Reflectivity -> jose value
    %inparam.sol_cD = 0; % Diffuse Reflectivity -> jose value is assumed to be 0 by heat program. Maybe value I cannot use
    
    verb = 0;

    % Initialize variables to store maximum values
    max_roll = 0;
    max_pitch = 0;
    max_yaw = 0;
    
    max_roll_unstable = 0; % roll
    max_pitch_unstable = 0; % pitch
    max_yaw_unstable = 0; % yaw

    max_p_roll = 0; % torques with propulsion
    max_p_pitch = 0;
    max_p_yaw = 0;
    
    max_p_roll_unstable = 0; % roll
    max_p_pitch_unstable = 0; % pitch
    max_p_yaw_unstable = 0; % yaw

    % average but over all positions!
    max_avg_roll = 0;
    max_avg_pitch = 0;
    max_avg_yaw = 0;

    n_fib = 200;
    all_roll_list = zeros(n_fib);
    all_pitch_list = zeros(n_fib);
    all_yaw_list = zeros(n_fib);

    max_all_roll = 0;
    max_all_pitch = 0;
    max_all_yaw = 0;


    % only has to be calculated once but for clarity its put here
    xg = sphere_fibonacci_grid_points(n_fib);

    % function to calculate stuff
    for cg_i = 1:numel(cg_s)
        cg = cell2mat(cg_s(cg_i));
        %display(cell2mat(cg))
        % Import model
        [modOut_imp] = ADBSatImport(modIn, modOut, verb, cg);
        
        % Nominal case
        for aos_j = -5:1:5  % Angle of sideslip [deg]
            for aoa_k = -1:0.1:1   % Angle of attack [deg] switched because blender wtf
                % Coefficient Calculation
                % fileOut = calc_coeff(modOut, resOut, deg2rad(aoa_deg), deg2rad(aos_deg), inparam, shadow, solar, 1, 0); 
                fileOut = calc_coeff(modOut_imp, resOut, deg2rad(aos_j), deg2rad(aoa_k), inparam, shadow, solar, 1, 0); 
                Cm_B = load(fileOut, "Cm_B").Cm_B;
                % Converter to aerodynamic torque!!
                rho = 5.33e-11;
                v_sat = 7729; % drago m/s
                M_b = 0.5* rho * Cm_B * load(fileOut, "AreaRef").AreaRef * v_sat^2;  % should be the actual torque calculation TODO

                % max_roll = max(max_roll, abs(M_b(1)));
                % max_pitch = max(max_pitch, abs(M_b(3))); % 3!!!!!!! around z here because either wrong model or wrong program
                % max_yaw = max(max_yaw, abs(M_b(2))); 
                % 
                % % check instability for pitch and yaw
                % if aoa_k*Cm_B(3) > 0
                %     max_pitch_unstable = max(max_pitch_unstable, abs(M_b(3)));
                % end
                % 
                % if aos_j*Cm_B(2) > 0
                %     max_yaw_unstable = max(max_yaw_unstable, abs(M_b(2)));
                % end

                % New setup

                [t_roll_torque, t_pitch_torque, t_yaw_torque, t_pitch_unstable, t_yaw_unstable, roll_torque, pitch_torque, yaw_torque, pitch_unstable, yaw_unstable] = torque_w_thrust(M_b, cg, aoa_k, aos_j);
                
                max_roll = max(max_roll, roll_torque);
                max_pitch = max(max_pitch, pitch_torque); % 3!!!!!!! around z here because either wrong model or wrong program
                max_yaw = max(max_yaw, yaw_torque); 
                max_pitch_unstable = max(max_pitch_unstable, pitch_unstable);
                max_yaw_unstable = max(max_yaw_unstable, yaw_unstable);
                % 
                max_p_roll = max(max_p_roll, t_roll_torque);
                max_p_pitch = max(max_p_pitch, t_pitch_torque); % 3!!!!!!! around z here because either wrong model or wrong program
                max_p_yaw = max(max_p_yaw, t_yaw_torque); 
                max_p_pitch_unstable = max(max_p_pitch_unstable, t_pitch_unstable);
                max_p_yaw_unstable = max(max_p_yaw_unstable, t_yaw_unstable);
            end
        end

        % Avg case
        for j = 1:n_fib
            [aos, aoa] = aos_aoa_converter(xg(j));
            fileOut = calc_coeff(modOut_imp, resOut, deg2rad(aos), deg2rad(aoa), inparam, shadow, solar, 1, 0); 
            Cm_B = load(fileOut, "Cm_B").Cm_B;
            % Converter to aerodynamic torque!!
            rho = 5.33e-11;
            v_sat = 7729; % drago m/s
            M_b = 0.5* rho * Cm_B * load(fileOut, "AreaRef").AreaRef * v_sat^2;  % should be the actual torque calculation TODO
            [roll_torque_all, pitch_torque_all, yaw_torque_all, ~, ~] = secondary_torques(M_b, aoa, aos); % dont really need unstable
            all_roll_list(j) = roll_torque_all;
            all_pitch_list(j) = pitch_torque_all;
            all_yaw_list(j) = yaw_torque_all;

            max_all_roll = max(max_all_roll, roll_torque_all);
            max_all_pitch = max(max_all_pitch, pitch_torque_all);
            max_all_yaw = max(max_all_yaw, yaw_torque_all);
        end
        max_avg_roll = max(max_avg_roll, mean(all_roll_list));
        max_avg_pitch = max(max_avg_pitch, mean(all_pitch_list));
        max_avg_yaw = max(max_avg_yaw, mean(all_yaw_list));
    end
end 


function [roll_torque, pitch_torque, yaw_torque, pitch_unstable, yaw_unstable] = secondary_torques(M_b, aoa, aos)
    pitch_torque = 0;
    yaw_torque = 0;
    % Gravity gradient torque (max 30 deg)
    Izz = 1.4;
    Iyy = 1.3;
    grav_constant = 3 / 2 * 3.986004418e14 / (6671800 ^ 3) ;
    %grav_torque = grav_constant*[(Izz - Iyy) * sin(2 * 30/180*pi), 0, 0]; % (Izz - Ixx) * sin(-2 * rot_angle[1])
    grav_torque_roll = grav_constant*(Izz - Iyy) * sin(2 * 30/180*pi);
    
    % Sun torque
    l_s_torque = 0.35;
    q_sun_max = 1414;
    a_albedo = 0.14;  % beta angle dependent 0.14 for < 30 0.19 above
    c_light = 299792458;  % [m/s] https://physics.nist.gov/cgi-bin/cuu/Value?c
    a_reflectivity = 0.9;
    a_frontal = 1.7;
    rad_torque = q_sun_max*(1+a_albedo)/c_light*a_frontal*(1+a_reflectivity)*l_s_torque;  % has to add in such a way that it makes it larger

    if aoa*M_b(3) > 0
        pitch_unstable = abs(M_b(3)) + rad_torque;
    else
        % create a unstable torque
        pitch_unstable = max(rad_torque - abs(M_b(3)), 0);
        % create a max torque
        pitch_torque = abs(M_b(3)) + rad_torque;
    end

    if aos*M_b(2) > 0
        yaw_unstable = abs(M_b(2)) + rad_torque;
    else
        yaw_unstable = max(rad_torque - abs(M_b(2)), 0);
        % create a max torque
        yaw_torque = abs(M_b(2)) + rad_torque;
    end
    
    sn_roll_t = sign(M_b(1));
    if sn_roll_t == 0  % when 0 torque can be added to any direction
        sn_roll_t  = 1;
    end
    roll_torque = abs(M_b(1) + sn_roll_t*(grav_torque_roll + rad_torque)); 
 
end




function [t_roll_torque, t_pitch_torque, t_yaw_torque, t_pitch_unstable, t_yaw_unstable, roll_torque, pitch_torque, yaw_torque, pitch_unstable, yaw_unstable] = torque_w_thrust(M_b, cg, aoa, aos)
    % Propulsion torque
    random_angle_offset = 0/180*pi;
    % minimises torque
    rad_offset = [9.63/180*pi, 0.1/180*pi];  % local angles compared to body angles, (roll is disregarded) first is a pitch angle and seconds is a yaw like angle

    prop_point = [0.18, 0, 0.146];  
    % if thruster in z direction is higher than cg then res pitch should be
    % positive (be aware that is third entry) (assuming that thrust is
    % aligned with flight direction
    % if in blender z direction moment with right hand rule is used that is 
    % actually negative
       
    % remember cg is y z x
    thrust = 8*10^-3;

    a = cg(3) - prop_point(1); % dk, i dont think it matters for now
    b = prop_point(2) - cg(1); % based on the text above
    c = prop_point(3) - cg(2);  % should be positive if z_point thruster is higher
    [roll_torque, pitch_torque, yaw_torque, pitch_unstable, yaw_unstable] = secondary_torques(M_b, aoa, aos);

    % add to M_b
    p_torque_roll = thrust*(b*sin(rad_offset(1)) - c*sin(rad_offset(2)));  % sin low
    p_torque_pitch = thrust*(c*cos(rad_offset(1))*cos(rad_offset(2)) - a*sin(rad_offset(1)));  % 
    %display(p_torque_pitch)
    p_torque_yaw = thrust*(a*sin(rad_offset(2)) - b*cos(rad_offset(1))*cos(rad_offset(2)));  % 

    M_b(1) = p_torque_roll + M_b(1);
    M_b(3) = p_torque_pitch + M_b(3); % yes 3
    M_b(2) = p_torque_yaw + M_b(2);

    [t_roll_torque, t_pitch_torque, t_yaw_torque, t_pitch_unstable, t_yaw_unstable] = secondary_torques(M_b, aoa, aos);
    t_roll_torque = t_roll_torque + abs(c*random_angle_offset);
    t_pitch_torque = t_pitch_torque + abs(a*random_angle_offset);
    t_pitch_unstable = t_pitch_unstable + abs(a*random_angle_offset);
    t_yaw_torque = t_yaw_torque + abs(a*random_angle_offset);
    t_yaw_unstable = t_yaw_unstable + abs(a*random_angle_offset);


    % p_torque_roll = b*rad_offset - c*rad_offset - c*random_angle_offset;  % 
    % p_torque_pitch = c*thrust -a*rad_offset - a*random_angle_offset;  % 
    % p_torque_yaw = a*rad_offset + a*random_angle_offset - b*thrust;  % 
end

function xg = sphere_fibonacci_grid_points ( ng )

%*****************************************************************************80
%
%% sphere_fibonacci_grid_points(): Fibonacci spiral gridpoints on a sphere.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    27 April 2015
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Richard Swinbank, James Purser,
%    Fibonacci grids: A novel approach to global modelling,
%    Quarterly Journal of the Royal Meteorological Society,
%    Volume 132, Number 619, July 2006 Part B, pages 1769-1793.
%
%  Input:
%
%    integer NG, the number of points.
%
%  Output:
%
%    real XG(N,3), the grid points.
%
  phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0;

  i = ( - ( ng - 1 ) : 2 : ( ng - 1 ) )';
  theta = 2 * pi * i / phi;
  sphi = i / ng;
  cphi = sqrt ( ( ng + i ) .* ( ng - i ) ) / ng;

  xg = zeros ( ng, 3 );

  xg(1:ng,1) = cphi .* sin ( theta );
  xg(1:ng,2) = cphi .* cos ( theta );
  xg(1:ng,3) = sphi;
   
end


function [aos, aoa] = aos_aoa_converter(x)
    % two circles
    % s_vector = [-1, 0, 0]  % xyz
    [azimuth, elevation, ~] = cart2sph(x(1), x(2), x(3));
    aos = azimuth;
    aoa = elevation;
end
% zero torque at precise cg 0.1765479 and angle 9.63194

%del = 0;

% aoa_deg = 5; % Angle of attack [deg] but not really this is the original definiton
% aos_deg = 1; % Angle of sideslip [deg] but not really
% % % Coefficient Calculation
% fileOut = calc_coeff(modOut, resOut, deg2rad(aoa_deg), deg2rad(aos_deg), inparam, shadow, solar, 1, 0); 
% 
% % Plot surface distribution
% if verb && ~del
%     plot_surfq(fileOut, modOut, aoa_deg(1), aos_deg(1), 'cp');
% end
%------------ END CODE -----------%