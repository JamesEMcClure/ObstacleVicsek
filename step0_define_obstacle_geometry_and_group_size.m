clearvars; close all; 
clc;

% This script is to define the geometry of the space you want your
% particles to live in. Define geometries with circular obstacles in a cell
% array, where each cell has the centers and radii of the obstacles in that
% case.  The script then generates some geometric descriptors of the space
% and saves them, plots the geometries for  you to see, and computes the
% number of agents you need for each case to keep the density fixed.

% Define the side length of square domain.
L = 10;      

% Define the particle density you want. 
rho = 2;  

%%%% Define the obstacle cases
obs_center{1} = [0,0]; obs_radii{1} = [0.1];    % case 1 
obs_center{2} = [0,0]; obs_radii{2} = [0.5];    % case 2 
obs_center{3} = [0,0]; obs_radii{3} = [1];      % case 3 
obs_center{4} = [0,0]; obs_radii{4} = [2];      % case 4 
obs_center{5} = [0,0]; obs_radii{5} = [4];      % case 5 

obs_center{6} = [-2,-2; 2,2]; obs_radii{6} = [0.1 0.1];     % case 6 
obs_center{7} = [-2,-2; 2,2]; obs_radii{7} = [0.5,0.5];     % case 7 
obs_center{8} = [-2,-2; 2,2]; obs_radii{8} = [1,1];         % case 8 
obs_center{9} = [-2,-2; 2,2]; obs_radii{9} = [1.5,1.5];     % case 9 
obs_center{10} = [-2,-2; 2,2]; obs_radii{10} = [2.5 2.5];   % case 10 

obs_center{11} = [-3,-3; 1.5,3; -1,0; -3,3; 1.5,-2; 3,-4]; 
obs_radii{11} = [0.1,0.1,0.1,0.1,0.1,0.1];                  % case 11
obs_center{12} = [-3,-3; 1.5,3; -1,0; -3,3; 1.5,-2; 3,-4]; 
obs_radii{12} = [0.25,0.25,0.25,0.25,0.25,0.25];            % case 12
obs_center{13} = [-3,-3; 1.5,3; -1,0; -3,3; 1.5,-2; 3,-4]; 
obs_radii{13} = [0.5,0.5,0.5,0.5,0.5,0.5];                  % case 13
obs_center{14} = [-3,-3; 1.5,3; -1,0; -3,3; 1.5,-2; 3,-4]; 
obs_radii{14} = [1,1,1,1,1,1];                              % case 14
obs_center{15} = [-3,-3; 1.5,3; -1,0; -3,3; 1.5,-2; 3,-4]; 
obs_radii{15} = [1.5,1.5,1.5,1.5,1,1];                      % case 15 

%%%% Compute area and perimeter of obstacle regions for each case

g = zeros(15,3);    % initialize data

for ii=1:15
    A = sum(pi*obs_radii{ii}.^2);
    P = sum(2*pi*obs_radii{ii});
    g(ii,:) = [A, P, (P^2)/A];
end


%%%% Compute the number of particles for each case to get density = rho

% N = density*(total area - obstacle area)
N = ceil(rho*(100-g(:,1)));

% Augment the N vector for two cases with no obstacles for when we generate data in the control conditions
N = [N; rho*L^2; rho*L^2];      


%%%% Plot a figure of all the cases

figure
for jj=1:5
    subplot(3,5,jj)
    for ii=1:length(obs_radii{jj})
        hold on, box on
        t = (0:1/100:1)*2*pi;
        x = obs_radii{jj}(ii)*cos(t)+obs_center{jj}(ii,1);
        y = obs_radii{jj}(ii)*sin(t)+obs_center{jj}(ii,2);
        fill(x,y,'k')
    end
    xlim([-L/2,L/2]),ylim([-L/2,L/2])
    axis square
    title(['Case ',num2str(jj)])
    
    kk=jj+5;
    subplot(3,5,kk)
    for ii=1:length(obs_radii{kk})
        hold on, box on
        t = (0:1/100:1)*2*pi;
        x = obs_radii{kk}(ii)*cos(t)+obs_center{kk}(ii,1);
        y = obs_radii{kk}(ii)*sin(t)+obs_center{kk}(ii,2);
        fill(x,y,'k')
    end
    xlim([-L/2,L/2]),ylim([-L/2,L/2])
    axis square
    title(['Case ',num2str(kk)])
    
    ll=jj+10;
    subplot(3,5,ll)
    for ii=1:length(obs_radii{11})
        hold on, box on
        t = (0:1/100:1)*2*pi;
        x = obs_radii{ll}(ii)*cos(t)+obs_center{ll}(ii,1);
        y = obs_radii{ll}(ii)*sin(t)+obs_center{ll}(ii,2);
        fill(x,y,'k')
    end
    xlim([-L/2,L/2]),ylim([-L/2,L/2])
    axis square
    title(['Case ',num2str(ll)])
end


%%%% Uncomment below to save the data and the figure

save('setup_data.mat','obs_center','obs_radii','N','g')
% set(gcf,'PaperPosition',[0,0,8,6]); print('-dpdf','obstacles.pdf')