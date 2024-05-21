%%%%%%%%%%%%%%% 
% Load necessary scripts for the solver and custom scripts for various functions
clear;
addpath(genpath("../00_essential_codes/"));
%%%%%%%%%%%%%%%

dg_globals = ReadYaml("input.yml");


if dg_globals.primary_spin == 0    
    SchwarzschildDriver
else

    %%% for future implementation
    KerrDriver
end