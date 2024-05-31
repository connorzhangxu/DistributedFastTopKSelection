%% Two time-scale distributed quantile estimatin
function [Error_Q]=DistributedQuantileEstimation_SGD(x,p,A,alpha0,beta0,tau1,tau2,n_iteration,Delta,loss)
D=diag(sum(A));
L=D-A;
N=length(x);
Error_Q=zeros(N,1);
[y,~]=sort(x,'ascend');
w=x;
i_iteration=1;
alpha=alpha0/(1+i_iteration)^tau1;
beta=beta0/(1+i_iteration)^tau2;
W=eye(N)-beta*L;
while i_iteration <=n_iteration %&& max_value> Delta/2

    dfw=(w>=x)-p;

    w=W*w-alpha*dfw;
    if strcmp(loss,'l2')
        Error_Q(i_iteration)=norm(w-y(floor(p*N)+1))^2/N;
    elseif strcmp(loss,'l1')
        Error_Q(i_iteration)=norm(w-y(floor(p*N)+1),1)/N;
    elseif strcmp(loss,'inf')
        Error_Q(i_iteration)=max(abs(w-y(floor(p*N)+1)));
    end
    i_iteration=i_iteration+1
end