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
inparam.Tw = 250; % Wall Temperature [K] -> jose

inparam.alpha = 1; % Accommodation (altitude dependent), 0.9 is more limiting from calculations
inc = 96.67; %deg
% 238 367 average for 65 65 

% for lat = -90:1:90 
%     for lon = -90:1:90
env = [alt*1e3, inc/2, 0, 106, 0, 238, 367, ones(1,7)*3, 0]; % Environment variables
% Environment Calculations
inparam = environment(inparam, env(1),env(2),env(3),env(4),env(5),env(6),env(7),env(8:14),env(15));

cg_s = {[0, 0, 0]};  % array of cg's y, z, x for blender its weird how it works

% end

modGeneral = 'Verification_Tests';
modNames = {'Block-Yaw', 'Block-Pitch', 'Block-Roll'};
get_config_torque_csv(cg_s, inparam, modNames, modGeneral)



function get_config_torque_csv(cg_s, inparam, modNames, modGeneral)
    ADBSat_path = ADBSat_dynpath;
    configs_torque_calculator(ADBSat_path, modNames, modGeneral, cg_s, inparam);
end

function configs_torque_calculator(ADBSat_path, modNames, modGeneral, cg_s, inparam)

    for modName_el = modNames
        modName = char(modName_el);
        display(modName)
        modIn = fullfile(ADBSat_path,'inou','obj_files','Configs', modGeneral,[modName,'.obj']);
        display(modIn)
        modOut = fullfile(ADBSat_path,'inou','models');
        resOut = fullfile(ADBSat_path,'inou','results', modName);
    
        calc_adcs_torques(cg_s, inparam, modIn, modOut, resOut);


    end
end

function calc_adcs_torques(cg_s, inparam, modIn, modOut, resOut)
    
    shadow = 1;

    solar = 0;
    %inparam.sol_cR = 0.77; % Specular Reflectivity -> jose value
    %inparam.sol_cD = 0; % Diffuse Reflectivity -> jose value is assumed to be 0 by heat program. Maybe value I cannot use
    
    verb = 0;

    % function to calculate stuff
    for cg_i = 1:numel(cg_s)
        cg = cell2mat(cg_s(cg_i));
        %display(cell2mat(cg))
        % Import model
        [modOut_imp] = ADBSatImport(modIn, modOut, verb, cg);
        aos_j = 0;
        aoa_k = 0;
        % Nominal case

        % Coefficient Calculation
        % fileOut = calc_coeff(modOut, resOut, deg2rad(aoa_deg), deg2rad(aos_deg), inparam, shadow, solar, 1, 0); 
        fileOut = calc_coeff(modOut_imp, resOut, deg2rad(aos_j), deg2rad(aoa_k), inparam, shadow, solar, 1, 0); 
        display(fileOut)
        Cm_B = load(fileOut, "Cm_B").Cm_B;
        % Converter to aerodynamic torque!!
        rho = 5.33e-11;
        v_sat = 7729; % drago m/s
        display(Cm_B, 'Cm_B')
        display(load(fileOut, "AreaRef").AreaRef, 'A_ref')
        M_b = 0.5* rho * Cm_B * load(fileOut, "AreaRef").AreaRef * v_sat^2 * load(fileOut, "LenRef").LenRef;  % should be the actual torque calculation TODO
        display(M_b(1), 'MB')  % roll
        display(M_b(3))  % pitch
        display(M_b(2))  % yaw

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

        %plot_surfq(fileOut, modOut, 0, 0, 'cp');

        % New setup

        % Avg case
    end
end 






%------------ END CODE -----------%