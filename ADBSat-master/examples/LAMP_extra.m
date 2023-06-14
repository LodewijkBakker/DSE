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

cg_s = {[0, 0.148, 0.002], [0, 0.158, 0.002], [0, 0.167, 0.002]};  % array of cg's y, z, x for blender its weird how it works

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


modNames = {'straight_combined_LAMP_0', 'straight_combined_LAMP_30', 'straight_combined_LAMP_60', 'straight_combined_LAMP_90'};
config_1_loc = 'C:\\Users\\Lodewijk\\Documents\\University_Subjects\\DSE\\ADBSat-master\\inou\\results\\config1.csv';

get_config_torque_csv(config_1_loc, cg_s, inparam, modNames)

function get_config_torque_csv(config_loc, cg_s, inparam, modNames)
    ADBSat_path = ADBSat_dynpath;
    config_f_id = fopen(config_loc,'w');
    fprintf(config_f_id, 'case, max torque roll, max torque pitch, max torque yaw, max unstable torque pitch, max unstable torque yaw\n');
    [abs_t_m_r, abs_t_m_p, abs_t_m_y, abs_t_m_p_u, abs_t_m_y_u] = configs_torque_calculator(ADBSat_path, modNames, cg_s, inparam, config_f_id);
    fprintf(config_f_id, 'absolute maximums\n');
    formatEnd = '%8.8f, %8.8f, %8.8f, %8.8f, %8.8f \n';
    fprintf(config_f_id, formatEnd, abs_t_m_r, abs_t_m_p, abs_t_m_y, abs_t_m_p_u, abs_t_m_y_u);
    fclose(config_f_id);
end

function [abs_t_m_r, abs_t_m_p, abs_t_m_y, abs_t_m_p_u, abs_t_m_y_u] = configs_torque_calculator(ADBSat_path, modNames, cg_s, inparam, config_1_f_id)
    abs_t_m_r = 0; % absolute torque max roll
    abs_t_m_p = 0; % absolute torque max pitch
    abs_t_m_y = 0; % absolute torque max yaw
    
    abs_t_m_p_u = 0; % absolute torque max pitch unstable
    abs_t_m_y_u = 0; % absolute torque max yaw unstable

    formatRes = '%s, %8.8f, %8.8f, %8.8f, %8.8f, %8.8f \n';

    for modName_el = modNames
        modName = char(modName_el);
        display(modName)
        modIn = fullfile(ADBSat_path,'inou','obj_files',[modName,'.obj']);
        modOut = fullfile(ADBSat_path,'inou','models');
        resOut = fullfile(ADBSat_path,'inou','results', modName);
    
        [t_m_r, t_m_p, t_m_y, t_m_p_u, t_m_y_u] = calc_adcs_torques(cg_s, inparam, modIn, modOut, resOut);
        % display(t_m_r)
        % display(t_m_p)
        % display(t_m_y)
        % display(t_m_p_u)
        % display(t_m_y_u)
    
        fprintf(config_1_f_id, formatRes, modName, t_m_r, t_m_p, t_m_y, t_m_p_u, t_m_y_u);

        abs_t_m_r = max(abs_t_m_r, t_m_r);
        abs_t_m_p = max(abs_t_m_p, t_m_p); % 3!!!!!!! around z here because either wrong model or wrong program
        abs_t_m_y = max(abs_t_m_y, t_m_y); 
            
        abs_t_m_p_u = max(abs_t_m_p_u, t_m_p_u);
        abs_t_m_y_u = max(abs_t_m_y_u, t_m_y_u);

    end
end

function [max_roll, max_pitch, max_yaw, max_pitch_unstable, max_yaw_unstable] = calc_adcs_torques(cg_s, inparam, modIn, modOut, resOut)
    
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

    % function to calculate stuff
    for cg_i = 1:numel(cg_s)
        cg = cell2mat(cg_s(cg_i));
        %display(cell2mat(cg))
        % Import model
        [modOut_imp] = ADBSatImport(modIn, modOut, verb, cg);
        
        for aos_j = -5:1:5  % Angle of sideslip [deg]
            for aoa_k = -1:0.1:1   % Angle of attack [deg] switched because blender wtf
                % Coefficient Calculation
                % fileOut = calc_coeff(modOut, resOut, deg2rad(aoa_deg), deg2rad(aos_deg), inparam, shadow, solar, 1, 0); 
                fileOut = calc_coeff(modOut_imp, resOut, deg2rad(aos_j), deg2rad(aoa_k), inparam, shadow, solar, 1, 0); 
                Cm_B = load(fileOut, "Cm_B").Cm_B;
                M_b = Cm_B * load(fileOut, "AreaRef").AreaRef;
                %display(Cm_B(1))
                max_roll = max(max_roll, abs(M_b(1)));
                max_pitch = max(max_pitch, abs(M_b(3))); % 3!!!!!!! around z here because either wrong model or wrong program
                max_yaw = max(max_yaw, abs(M_b(2))); 
        
                % check instability for pitch and yaw
                if aoa_k*Cm_B(3) > 0
                    max_pitch_unstable = max(max_pitch_unstable, abs(M_b(3)));
                end
        
                if aos_j*Cm_B(2) > 0
                    max_yaw_unstable = max(max_yaw_unstable, abs(M_b(2)));
                end
            end
        end
    end
end 


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