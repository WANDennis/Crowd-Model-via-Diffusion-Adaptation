function [Disagree,MSE,MSD]=nocoop(R,r,times)
warn = warning ('off','all');

% range, startpoint and target
Xmin=-40;           %range of image
Ymin=-10;
Xmax=40;
Ymax=10;
start_point=[-40 ; 0]; % start point of all nodes
var=5;% point generate variance
wt=[200 ; 10]; % target point location

% coefficients
N=50;      % no. of nodes
dt=0.5;     % time increment / time step
% times=100; %Number of iterations
% R=5;   % neighborhood radius
% r=3; % node distance limit
r_min=3;% minimum tolerant distance
d=1; % distance of population from the walls
d_high=3; % distance of population from the walls (green)
chord_threshold=34; % wall width (chord length) change threshold 
x_chord_begin=20;
x_chord_end=140;
tail_l=0.5;
dmax=1.2;

%velocity parameters
lamda=0.5;  
alpha=2;
alpha_max=4;
beta=2;
gamma=0.5; 
mu=0.5; % step size
kappa=0.01;

% define wall models
syms P1(z) P2(z)
P1(z)= 0.008*(z)^2 +100;
P2(z)= -0.008*(z)^2 +-100;
Pd1=diff(P1);% derivative of wall equations
Pd2=diff(P2);
steps1=200;
steps2=steps1;
z1=linspace(Xmin,Xmax,steps1).';
z2=linspace(Xmin,Xmax,steps2).';
yP1=P1(z1);
yP2=P2(z2);
curvexy1=[z1 yP1];
curvexy2=[z2 yP2];
syms y1 x1 y2 x2 cx cy real

% boundaries of wall (red)
syms P1_limit P2_limit
P1_limit(z)= P1(z) -d;
P2_limit(z)= P2(z) +d;
yP1_limit=P1_limit(z1);
yP2_limit=P2_limit(z2);
curvexy_limit1=[z1 yP1_limit];
curvexy_limit2=[z2 yP2_limit];

% boundaries of wall (green)
syms P1_limit_green P2_limit_green
P1_limit_green(z)= P1(z) -d_high;
P2_limit_green(z)= P2(z) +d_high;
yP1_limit_green=P1_limit_green(z1);
yP2_limit_green=P2_limit_green(z2);
curvexy_limit3=[z1 yP1_limit_green];
curvexy_limit4=[z2 yP2_limit_green];

% directing 
interval=30;
steps_direct=linspace(start_point(1),wt(1),interval).';% x
directx=[steps_direct(2:end-1)];
midpoint=double((P1(directx)+P2(directx))/2); % y
directxy=[directx.'; midpoint.'];
directxy=[directxy,wt];
directxy=wt;
directx=directxy(1,:);
directy=directxy(2,:);

% generating the initial population (distributed in acircle with radius 'var')
x=zeros(2,N,times+1); % location of all nodes (uniformly distribute in a circle with radius var)
rng(7);
radius = var*sqrt(rand(1,N)); % Using square root here ensures distribution uniformity by area
rng(8);
t = 2*pi*rand(1,N);
x_cor = radius.*cos(t);
y_cor = radius.*sin(t);
x(:,:,1)=[start_point(1)+x_cor;start_point(2)+y_cor];
x(:,:,2)=[start_point(1)+x_cor;start_point(2)+y_cor];

% define all variables
u_nocoop=zeros(2,N,times+1); % velocity of all nodes
ua=zeros(2,N,times+1); % heading target velocity of all nodes
ub=zeros(2,N,times+1);% reunion velocity
uc=zeros(2,N,times+1); % coherent motion velocity
delta=zeros(2,N,times+1);
ug_nocoop=zeros(2,N,times); % gravity center velocity
ug_ture_nocoop=zeros(2,times);                                   %%%%%%%%%%%%%%%%%%%%%
dist=zeros(N,times); % distance from each node to estimated target
dist_ture=zeros(N,times); % distance from each node to estimated target %%%%%%%%%%%%%
phi=zeros(2,N,times+1);
varphi=zeros(2,N,times+1);   %%%%%%%%%%%%%%%%%%%%%
flag_edge=zeros(N,times);
chord_length=zeros(N,times);
chord_centre=zeros(2,N,times);
r_change=r*ones(N,times);
w_nocoop=zeros(2,N,times);   % estimated target locatrion for every k
alpha_change=alpha*ones(N,times+1);
vind_all_include=cell(N);
tailxy=zeros(2,N,times+1);
q=zeros(2,N,times);  %%%%%%%%%%%%%%%%%%
w_true=zeros(2,times);  %%%%%%%%%%%%%%
noise_variance=zeros(N,times+1);%%%%%%%%%%%
xg=zeros(2,times);
dist_g=zeros(N,times);

% initialize velocity 
rng(6);
randv=randn(2,N);
u_nocoop(:,:,2)=randv./vecnorm(randv);

% at time t, location is updated by u(:,k,t+1)
for t=2:times

    w_true(:,t)=directxy;   %%%%%%%%%%%%%%%%%%%%%%

    dist_nodek=zeros(N,N); % distance between every two nodes (l,k)
    
    % center of gravity
    xg(:,t)=mean(x(:,:,t),2);
    % node cloest to center of gravity
    dist_xg=vecnorm(x(:,:,t)-xg(:,t)).';
    [~, vind_g]=min(dist_xg);

    for k=1:N

        % update neighbourhood
        dist_nodek(:,k)=vecnorm(x(:,:,t)-x(:,k,t)).'; % distance between k and its neighbour
        vind=find(dist_nodek(:,k)<=R & dist_nodek(:,k)>0); % index of neighbourhood (doesn't include itself)
        vind_all_include{k}=find(dist_nodek(:,k)<=R); % index of  neighbourhood (include itself)
        vind_nonneigh=find(dist_nodek(:,k)>R);% index of outside neighbourhood 
        Nk=length(vind_all_include{k}); % number of neighbours of k
        
        % store vin_g data
        if k == vind_g
            selected=dist_nodek(dist_nodek(:,k)<=R,k);
%             [~,I]=sort(dist_nodek(:,k),'ascend');
%             vind_g_neigh=I(2:11);
            dist_g(1:length(selected),t)=selected;
        end



        % calculate distance from nodes to estimated target
        % update w (estimated w for every k)   %%%%%%%%%%%%%%%%%%%
        w_nocoop(:,k,t)=varphi(:,k,t);
%       w(:,k,t)=w_true(:,t); %%%% set default
        % dist between each node and target
        dist(:,t)=vecnorm(w_nocoop(:,:,t)-x(:,:,t)).'; 

        ua(:,k,t+1)=(w_nocoop(:,k,t)-x(:,k,t))./dist(k,t);
        

        % update gravity center velocity
        ug_nocoop(:,k,t)=phi(:,k,t);
        % update true gravity center velocity   %%%%%%%%%%%%%%%%%%%%%%
        ug_ture_nocoop(:,t)=mean(u_nocoop(:,:,t),2);


%         % reunion 
%         % decide which edge
%         W=[u_nocoop(1,k,t)/norm(u_nocoop(:,k,t)) -u_nocoop(2,k,t)/norm(u_nocoop(:,k,t)); u_nocoop(2,k,t)/norm(u_nocoop(:,k,t)) u_nocoop(1,k,t)/norm(u_nocoop(1,k,t))]; 
%         xlk=zeros(2,N);
%         for l=vind  % just care about neighbor fish
%             xlk(:,l)=W'*(x(:,l,t)-x(:,k,t)); % (W times every relative location of fish k)
%         end
%         if xlk(1,:)==-abs(xlk(1,:)) % (all negative, no node l lies in the front of k) k is at front
%            flag_edge(k,t)=1;
%         elseif xlk(2,:)==abs(xlk(2,:)) % (all positive, all nodes l lie in left of node k) l(s) are all at left, k is at right edge
%            flag_edge(k,t)=2;
%         elseif xlk(2,:)==-abs(xlk(2,:)) % l(s) are all at right, k is at left edge
%            flag_edge(k,t)=3;
%         else 
%            flag_edge(k,t)=0;
%         end
       
        % update coherent motion velocity (avoid collisions)
        if Nk > 1
                delta(:,k,t)=(1/(Nk-1)) * sum((1-r_change(k,t)./dist_nodek(vind,k)).' .*  (x(:,vind,t)-x(:,k,t)),2); 
        end
        uc(:,k,t+1)=(1-lamda)*ug_nocoop(:,k,t)+gamma*delta(:,k,t);
        
        % update velocity of each node 
        u_nocoop(:,k,t+1)=lamda*(alpha_change(k,t+1)*ua(:,k,t+1)+beta*ub(:,k,t))+uc(:,k,t+1);

        % update phi (local center of gravity velocity elements)
        phi(:,k,t+1)=(1-mu)*ug_nocoop(:,k,t)+mu*u_nocoop(:,k,t+1);

        % update location
        x(:,k,t+1)=x(:,k,t)+dt*u_nocoop(:,k,t+1);
       
       
        w_true(:,t+1)=directxy;   %%%%%%%%%%%%%%%%%%%%%%
        dist_ture(:,t+1)=vecnorm(w_true(:,t+1)-x(:,:,t+1)).'; % dist between each node and target   %%%%%%%%%%%%%%%%%%%%
        % add noise  %%%%%%%%%%%%
        noise_variance(k,t+1)=sqrt(kappa*dist_ture(k,t+1).^2);
        rng('shuffle');
        q(:,k,t+1)= w_true(:,t+1)+noise_variance(k,t+1)*randn(2,1);
        % update varphi (target location elements) %%%%%%%%%%%%%%%%%
        varphi(:,k,t+1)= (1-mu)*w_nocoop(:,k,t) + mu*q(:,k,t+1);
    end


%     % set reunion velocity
%     % find the nearest edge (different from itself)
%     for k=1:N
%         flag_out_neigh=flag_edge(:,t);
%         flag_out_neigh(vind_all_include{k})=0;% flag for nodes outside the neighbour of node k
%         if flag_edge(k,t)==1 % front
%            tb =find(flag_out_neigh == 1); % the index of same edges
%            if length(tb) == 0
%                ub(:,k,t+1)=zeros(2,1);
%            else
%                dist_kl=dist_nodek(tb,k);    
%                [~, ind]= min(dist_kl);
%                l=tb(ind);
%                ub(:,k,t+1)= (x(:,l,t) - x(:,k,t))/norm(x(:,l,t) - x(:,k,t));
%            end
%        
%         elseif flag_edge(k,t)==2 % right
%            tb =find(flag_out_neigh == 2); % the index of same edges
%            if length(tb) == 0
%                ub(:,k,t+1)=zeros(2,1);
%            else
%                dist_kl=dist_nodek(tb,k);    
%                [~, ind]= min(dist_kl);
%                l=tb(ind);
%                ub(:,k,t+1)= (x(:,l,t) - x(:,k,t))/norm(x(:,l,t) - x(:,k,t));
%            end
% 
%         elseif flag_edge(k,t)==3 % left
%            tb =find(flag_out_neigh == 3); % the index of same edges
%            if length(tb) == 0
%                ub(:,k,t+1)=zeros(2,1);
%            else
%                dist_kl=dist_nodek(tb,k);    
%                [~, ind]= min(dist_kl);
%                l=tb(ind);
%                ub(:,k,t+1)= (x(:,l,t) - x(:,k,t))/norm(x(:,l,t) - x(:,k,t));
%            end
%         else 
%             ub(:,k,t+1)=zeros(2,1);
%         end
%     end 

end

Disagree=zeros(1,times-1);
MSE=zeros(1,times-1);
MSD=zeros(1,times-1);
for t=2:times
    Disagree(t-1)=mean(vecnorm(u_nocoop(:,:,t)-ug_ture_nocoop(:,t)).^2);
    MSE(t-1)=mean(vecnorm(ug_nocoop(:,:,t)-ug_ture_nocoop(:,t)).^2);
    MSD(t-1)=mean(vecnorm(w_nocoop(:,:,t)-w_true(:,t)).^2);
end