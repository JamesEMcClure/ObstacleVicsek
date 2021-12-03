clearvars; close all; clc;
set(0,'DefaultFigureWindowStyle','docked');

% This script is to generate data from the Vicsek model with periodic BCs and
% when obstacles are in the domain. Collisions with the obstacles are 
% reflective. It also generates data for two control conditions: 1) where 
% there are no obstacles and 2) where agents do not interact and only 
% respond to additive noise (i.e. random walkers).

%%%%%%%%%%%%%% Define parameters 
maxT = 1000;            % sim length
L = 10;                 % side length, must be >= 10 so obstacle choices below make sense
s = 0.15;               % particle speed
reps = 20;              % number of monte carlo replicates
ncases = 17;            % number of obstacle cases + two control cases

% alignment noise, range [0,1]
% etaV = [0.01 0.02 0.04 0.05 0.08 0.1 0.16 0.2 0.22 0.25 0.28 0.3 0.32 0.35 0.38 0.4 0.41 0.42 0.43 0.44 0.45 0.5];
etaV = [0.02];

% load number of particles and obstacle geometry
load('setup_data.mat');

% add path to find the Vicsek scripts
addpath('.\Vicsek_scripts');

%%%%%%%%%%%%%% Generate data
for irep = 1:reps
    
    %%%%%%%%%%%%%% Initialize data structures
    pos = cell(ncases);
    pos_nm = cell(ncases);
    vel = cell(ncases);
    
    for ieta = 1:length(etaV)
        eta = etaV(ieta);
        
        % run cases with obstacles
        for icase = [1 11]
            [auxpos,auxpos_nm,auxvel] = vicsek_periodic_obstacles(N(icase),L,maxT,s,eta,obs_center{icase},obs_radii{icase});
            pos{icase}(:,:,ieta) = auxpos;
            pos_nm{icase}(:,:,ieta) = auxpos_nm;
            vel{icase}(:,:,ieta) = auxvel;
        end
%         
%         % run control condition with no obstacles
%         [auxpos,auxpos_nm,auxvel] = vicsek_periodic(N(16),L,maxT,s,eta);
%         pos{16}(:,:,ieta) = auxpos;
%         pos_nm{16}(:,:,ieta) = auxpos_nm;
%         vel{16}(:,:,ieta) = auxvel;
%         
%         % run control condition with no interaction
%         [auxpos,auxpos_nm,auxvel] = vicsek_periodic_control(N(17),L,maxT,s,eta);
%         pos{17}(:,:,ieta) = auxpos;
%         pos_nm{17}(:,:,ieta) = auxpos_nm;
%         vel{17}(:,:,ieta) = auxvel;
%     
    end
%     irep
    %%%% Uncomment below to save the data 
%     save(['data/data_rep',num2str(irep),'.mat'],'pos','pos_nm','vel','maxT','N','L','s','etaV','reps');
end
