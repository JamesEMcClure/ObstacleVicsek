function [datap,datap_nm,datav] = vicsek_periodic_control(N,L,maxtime,s,eta)

% 2D Vicsek model with periodic BCs and NO INTERACTION

rZI = 1;            % radius of the zone of interaction
tau = 1.0;          % time step

% initialize to random positions and velocities
cx = L*(rand([N,1])-0.5);
cy = L*(rand([N,1])-0.5);
theta = 2*pi*rand([N,1]);
thetaHat=theta;
vx = s*cos(theta);
vy = s*sin(theta);

cx_nm = cx; cy_nm = cy;

% preallocate space for data
datap=zeros(maxtime,2*N);
datav=zeros(maxtime,2*N);

for itime = 1 : maxtime
    
    % initialize for thie time step
    cx1=[]; cy1=[]; vx1=[]; vy1=[];
    cxA=cx; cyA=cy;vxA=vx; vyA=vy;
    
    % find desired orientations 
    for ii=1:N
        noise = rand-0.5;
        thetaHat(ii) = atan2(vyA(ii),vxA(ii)) + eta*2*pi*noise;
    end
   
    % provisional update velocity and position
    vx1 = s*cos(thetaHat);
    vy1 = s*sin(thetaHat);
  
    cx1 = cx + vx1*tau;
    cy1 = cy + vy1*tau;
    
    cx1 = mod(cx1+L/2,L) - L/2;
    cy1 = mod(cy1+L/2,L) - L/2;
    
    cx_nm = cx_nm + vx1*tau; 
    cy_nm = cy_nm + vy1*tau;
    
    % real update
    cx=cx1; cy=cy1; vx=vx1; vy=vy1;
    
    % accumulate data
    dd=[cx,cy]'; dd=dd(:)'; ee=[vx, vy]'; ee=ee(:)';
    datap(itime,:)=dd; datav(itime,:)=ee;
        
    datap_nm(itime,:) = [cx_nm' cy_nm'];
end
end