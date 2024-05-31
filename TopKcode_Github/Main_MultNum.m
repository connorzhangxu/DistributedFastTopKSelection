%% Distributed smoothed quantile inference for multiple data points within each node

clc
clear
close all
addpath(genpath('./utils/'));
addpath(genpath('./algorithms/'));
addpath(genpath('./data/'));

N=1e5; % Number of data
n=1e3; % Number of nodes          
Delta = 0.1; %  resolution 
k=round(N*0.1);

n_iteration=3e4; % Maximum iteration

%% set random seed
seed=100;
rng(seed);

%% generate signal with resolution delta

x=round(randn(N,1)*sqrt(10)/Delta)*Delta; % data

data=RandAssignData(N,n,x);

%% generate Erdo Renyi random graph
NumEdges=3*n;
[A]=RandomGraphGeneration(n,NumEdges);
% figure
% plot(graph(A))
dmax=max(sum(A));
D=diag(sum(A));
L=D-A;
lambda=svd(L);

%% Main program
tau1=0;
tau2=0;
alpha0=0.003*Delta;
beta0=2/(lambda(1)+lambda(n-1));

p=(N-k+0.5)/N;
[y,~]=sort(x,'descend');
m_over=k-sum(x>y(k));
m_under=N-k-sum(x<y(k));
gm=min(m_over-0.5,m_under+0.5);
W=eye(n)-beta0*L;

% loss='l2';
% loss='l1';
loss='inf';
% smooth='Nesterov';
smooth='Convolution';
threshold=y(k);

Error_Q1=DistributedQuantileEstimation_SGD_MultNum(data,threshold,p,A,alpha0,beta0,tau1,tau2,n_iteration,Delta,loss);

h=Delta*0.1;
Error_Q2=DistributedQuantileEstimation_EXTRA_MultNum(data,threshold,p,A,beta0,h,n_iteration,Delta,loss,smooth);

%% Plot
figure
X=round(logspace(0,4.47,50));
loglog(X,Error_Q1(X),'-^','linewidth',2)
% loglog(Error_Q1,'-^','linewidth',2)
hold on
loglog(X,Error_Q2(X),'-s','linewidth',2)
% loglog(Error_Q2,'-s','linewidth',2)
loglog([1:n_iteration],Delta/2*ones(n_iteration,1),'k-','linewidth',1.5)

legend('DGD','EXTRA','LineWidth',1.5)
xlabel('$t$','interpreter','latex')
axis([0 3e4 1e-3 1.5e1])
if strcmp(loss,'l2')
    ylabel('$\|\mathbf{w}^t - \theta_k \mathbf{1}\|_{2}^2/N$','interpreter','latex')
elseif strcmp(loss,'max')
    ylabel('$\|\mathbf{w}^t - \theta_k \mathbf{1}\|_{\infty}$','interpreter','latex')
elseif strcmp(loss,'l1')
    ylabel('$\|\mathbf{w}^t - \theta_k \mathbf{1}\|_{1}/N$','interpreter','latex')
end

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6 4.5]);
set(gca,'FontName','times new roman','FontSize',16,'Layer','top','LineWidth',2);
