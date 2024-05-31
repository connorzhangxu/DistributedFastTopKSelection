%% Two time-scale distributed quantile estimatin
function [Error_Q]=DistributedQuantileEstimation_SGD_MultNum(data,threshold,p,A,alpha0,beta0,tau1,tau2,n_iteration,Delta,loss)
D=diag(sum(A));
L=D-A;
n=length(data);
Error_Q=zeros(n_iteration,1);
w=zeros(n,1);
for i=1:n
    w(i)=mean(data{i});
end
i_iteration=1;
while i_iteration <=n_iteration %&& max_value> Delta/2
    alpha=alpha0/(1+i_iteration)^tau1;
    beta=beta0/(1+i_iteration)^tau2;
    dfw=zeros(n,1);
    for i=1:n
        dfwt=(w(i)>=data{i})-p;
        dfw(i)=sum(dfwt);
    end

    W=eye(n)-beta*L;
    w=W*w-alpha*dfw;
    if strcmp(loss,'l2')
        Error_Q(i_iteration)=norm(w-threshold)^2/n;
    elseif strcmp(loss,'l1')
        Error_Q(i_iteration)=norm(w-threshold,1)/n;
    elseif strcmp(loss,'inf')
        Error_Q(i_iteration)=max(abs(w-threshold));
    end
    i_iteration=i_iteration+1
end