function [] = particle_filter_localization()
%PARTICLE_FILTER_LOCALIZATION Summary of this function goes here
%   Detailed explanation goes here

% -------------------------------------------------------------------------
% TASK for particle filter localization
% for robotic class in 2018 of ZJU

% Preparartion: 
% 1. you need to know how to code and debug in matlab
% 2. understand the theory of Monte Carlo

% Then complete the code by YOURSELF!
% -------------------------------------------------------------------------

close all;
clear all;

disp('Particle Filter program start!!')

%% initialization

time = 0;
endTime = 60; % second
global dt;
dt = 0.1; % second
 
nSteps = ceil((endTime - time)/dt);
  
localizer.time = [];
localizer.xEst = [];
localizer.xGnd = [];
localizer.xOdom = [];
localizer.z = [];
localizer.PEst=[];
localizer.u=[];

% Estimated State [x y yaw]'ä¼°è?¡å??
xEst=[0 0 0]';
% GroundTruth State????
xGnd = xEst;
% Odometry-only = Dead Reckoning  
xOdom = xGnd;

% Covariance Matrix for predict
Q=diag([0.1 0.1 toRadian(3)]).^2;
Q
% Covariance Matrix for observation
R=diag([1]).^2;% range:meter
R
% Simulation parameter??
global Qsigma
Qsigma=diag([0.1 toRadian(5)]).^2;
global Rsigma
Rsigma=diag([0.1]).^2;

% landmark position
landMarks=[10 0; 10 10; 0 15; -5 20];

% longest observation confined
MAX_RANGE=20;
% Num of particles, initialized
NP=50;
% Used in Resampling Step, a threshold
NTh=NP/2.0;

% particles produced
px=repmat(xEst,1,NP);
% weights of particles produced
pw=zeros(1,NP)+1/NP;


%% Main Loop 

for i=1 : nSteps
    
    time = time + dt;
    u=doControl(time);
    
    % do observation
    [z,xGnd,xOdom,u]=doObservation(xGnd, xOdom, u, landMarks, MAX_RANGE);
    
    for ip=1:NP
        
        % process every particle
        x=px(:,ip);
        w=pw(ip);
        
        % do motion model and random sampling
        x=doMotion(x, u)+sqrt(Q)*randn(3,1);
    
         % calculate inportance weight
        for iz=1:length(z(:,1))
            pz=norm(x(1:2)'-z(iz,2:3));
            dz=pz-z(iz,1);
            w=w*Gaussian(dz,0,sqrt(R));
        end
        px(:,ip)=x;
        pw(ip)=w;
        
    end
    
    pw=Normalization(pw,NP);
    [px,pw]=ResamplingStep(px,pw,NTh,NP);
    xEst=px*pw';
    
    % Simulation Result
    localizer.time=[localizer.time; time];
    localizer.xGnd=[localizer.xGnd; xGnd'];
    localizer.xOdom=[localizer.xOdom; xOdom'];
    localizer.xEst=[localizer.xEst;xEst'];
    localizer.u=[localizer.u; u'];
    
    % Animation (remove some flames)
    if rem(i,10)==0 
        hold off;
        arrow=0.5;
        for ip=1:NP
            quiver(px(1,ip),px(2,ip),arrow*cos(px(3,ip)),arrow*sin(px(3,ip)),'ok');hold on;
        end
        plot(localizer.xGnd(:,1),localizer.xGnd(:,2),'.b');hold on;
        plot(landMarks(:,1),landMarks(:,2),'pk','MarkerSize',10);hold on;
        if~isempty(z)
            for iz=1:length(z(:,1))
                ray=[xGnd(1:2)';z(iz,2:3)];
                plot(ray(:,1),ray(:,2),'-r');hold on;
            end
        end
        plot(localizer.xOdom(:,1),localizer.xOdom(:,2),'.k');hold on;
        plot(localizer.xEst(:,1),localizer.xEst(:,2),'.r');hold on;
        axis equal;
        grid on;
        drawnow;
    end
    
end

% draw the final results of localizer, compared to odometry & ground truth
drawResults(localizer);


end









%% Other functions

% degree to radian
function radian = toRadian(degree)
    radian = degree/180*pi;
end

function []=drawResults(localizer)
%Plot Result
 
    figure(1);
    hold off;
    x=[ localizer.xGnd(:,1:2) localizer.xEst(:,1:2)];
    set(gca, 'fontsize', 12, 'fontname', 'times');
    plot(x(:,1), x(:,2),'-.b','linewidth', 4); hold on;
    plot(x(:,3), x(:,4),'r','linewidth', 4); hold on;
    plot(localizer.xOdom(:,1), localizer.xOdom(:,2),'--k','linewidth', 4); hold on;

    title('Localization Result', 'fontsize', 12, 'fontname', 'times');
    xlabel('X (m)', 'fontsize', 12, 'fontname', 'times');
    ylabel('Y (m)', 'fontsize', 12, 'fontname', 'times');
    legend('Ground Truth','Particle Filter','Odometry Only');
    grid on;
    axis equal;

end

function [ u ] = doControl( time )
%DOCONTROL Summary of this function goes here
%   Detailed explanation goes here

    %Calc Input Parameter
    T=10; % [sec]

    % [V yawrate]
    V=1.0; % [m/s]
    yawrate = 5; % [deg/s]

    u =[ V*(1-exp(-time/T)) toRadian(yawrate)*(1-exp(-time/T))]';


end


%%  you need to complete

%% do Observation model 
function [z, xGnd, xOdom, u] = doObservation(xGnd, xOdom, u, landMarks, MAX_RANGE)
    global Qsigma;
    global Rsigma;
    
    % Gnd Truth and Odometry
    xGnd=doMotion(xGnd, u);% Ground Truth
    u=u+sqrt(Qsigma)*randn(2,1); % add noise randomly
    xOdom=doMotion(xOdom, u); % odometry only
    
    %Simulate Observation
    z=[];
    for iz=1:length(landMarks(:,1))
        d = norm(xGnd(1:2)' - landMarks(iz,:));%d = norm( -landMarks(:,1) )

        if d<MAX_RANGE 
            z=[z;[d+sqrt(Rsigma)*randn(1,1) landMarks(iz,:)]];   % add observation noise randomly
        end
    end
end


%% do Motion Model
function [ x ] = doMotion( x, u)
    global dt;

    Delta = [ dt*cos(x(3)) 0
              dt*sin(x(3)) 0
              0 dt];
    Identity = eye(3)
    x = Identity*x + Delta*u
      
end

%% Gauss function
function g = Gaussian(x,u,sigma)
    g= gaussmf(x,[sigma u])
end

%% Normalization 
function pw=Normalization(pw,NP)
    pwsum = sum(pw)
    for i = 1 : NP
        pw(i) = pw(i) / pwsum 
    end
end

%% Resampling
function [px,pw]=ResamplingStep(px,pw,NTh,NP)
    N_eff= 1/(pw*pw');
    if N_eff < NTh
        ppx=px
    for i=1:NP
      u=rand;
      pw_sum=0;
      for j=1:NP
          pw_sum = pw_sum + pw(j);
          if pw_sum >= u
              px(:,i) = ppx(:,j);
              break;
          end
      end
    end
    pw = zeros(1,NP)+1/NP;
    end
end

