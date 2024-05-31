%% Distributed smoothed quantile inference for large scale case

clc
clear
close all
addpath(genpath('./utils/'));
addpath(genpath('./algorithms/'));
addpath(genpath('./data/'));

N=1e4; % Number of nodes
n_iteration=2e4; % Maximum iteration
k=round(N*0.3);
%% set random seed
seed=10;
rng(seed);

%% generate signal with resolution delta
Delta = 0.1;
x=round(randn(N,1)*sqrt(10)/Delta)*Delta;

%% generate Erdo Renyi random graph
% NumEdges=5*N;
% [A]=RandomGraphGeneration(N,NumEdges);
% figure
% plot(graph(A))
% diameter
% dmax=max(sum(A));
% D=diag(sum(A));
% L=D-A;
% lambda=svd(L);

load('graph1e4.mat')
%% Main program
tau1=0;
tau2=0;
alpha0=0.04*Delta;
beta0=2/(lambda(1)+lambda(N-1));
p=(N-k+0.5)/N;
[y,~]=sort(x,'descend');
m_over=k-sum(x>y(k));
m_under=N-k-sum(x<y(k));
gm=min(m_over-0.5,m_under+0.5);
% W=eye(N)-beta0*L;
% Sigma=svd(W);

% loss='l2';
% loss='l1';
loss='inf';
% smooth='Nesterov';
smooth='Convolution';

Error_Q1=DistributedQuantileEstimation_SGD(x,p,A,alpha0,beta0,tau1,tau2,n_iteration,Delta,loss);
h=Delta*5;
Error_Q2=DistributedQuantileEstimation_EXTRA(x,p,A,beta0,h,n_iteration,Delta,loss,smooth);

%% Plot
figure
% X=round(logspace(0,4,60));
X=round(logspace(0,log10(n_iteration),60));
loglog(X,Error_Q1(X),'-^','linewidth',2)
% loglog(Error_Q1,'-s','linewidth',2)
hold on
loglog(X,Error_Q2(X),'-^','linewidth',2)
% loglog(Error_Q2,'-s','linewidth',2)
hold on
loglog([1:n_iteration],Delta/2*ones(n_iteration,1),'k-','linewidth',1.5)
axis([0 n_iteration 8e-4 1e1])

legend('DGD','EXTRA','LineWidth',1.5)
xlabel('$t$','interpreter','latex')
ylabel('$||\mathbf{w}^t-\theta_k \mathbf{1}||_\infty$','interpreter','latex')
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
