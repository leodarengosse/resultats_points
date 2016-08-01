function [] = run_plot_antex(ant_type, freq, save)
%% function [] = run_plot_antex(ant_type, freq, save)
%% 
%% Run plot antex
%% 
%% Clement Fontaine - 2014-01-13
%%
%% Input : 
%% - ant_type : antenna name, (IGS name)
%% - freq : frequency among ['G01','G02','R01','R02']
%% - savename : save name (optional)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% antex loading
[ATX_header, ATX_data] = load_antex('igs08.atx');

% plot antex
if nargin == 4
	savename = save;
else
	savename = '';
end

plot_antex(ATX_header, ATX_data, ant_type, freq, savename);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
