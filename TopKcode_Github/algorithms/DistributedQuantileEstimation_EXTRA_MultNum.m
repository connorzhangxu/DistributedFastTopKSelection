%% EXTRA
function [Error_Q]=DistributedQuantileEstimation_EXTRA_MultNum(data,threshold,p,A,beta0,h,n_iteration,Delta,loss,smooth_type)

D=diag(sum(A));
L=D-A;
n=length(data);
Error_Q=zeros(n_iteration,1);
v=zeros(n,1);
W=eye(n)-beta0*L;

w=zeros(n,1);
for i=1:n
    w(i)=mean(data{i});
end
% beta=1/h/sqrt(1-Sigma(2));
beta=1/h;
alpha=1/beta;
max_value=inf;
i_iteration=1;

while i_iteration <=n_iteration %&& max_value> Delta/2
    dfw=zeros(n,1);
    if strcmp(smooth_type,'Nesterov')
        % Nesterov's smoothing
        for i=1:n
            dfwt=((data{i}-w(i))<=h*(p-1))*(1-p)-((data{i}-w(i))>h*p)*p+(h*(p-1)<(data{i}-w(i)) & (data{i}-w(i))<=h*p).*(w(i)-data{i})/h;
            dfw(i)=sum(dfwt);
        end
    elseif strcmp(smooth_type,'Convolution')
        % Convolution smoothing
        for i=1:n
            dfwt=((data{i}-w(i))<=-h)*(1-p)-((data{i}-w(i))>h)*p+(-h<(data{i}-w(i)) & (data{i}-w(i))<=h).*((w(i)-data{i}+h)/2/h-p);
            dfw(i)=sum(dfwt);
        end
    else
        fprintf('Please input right smoothing type! Options include Nesterov and Convolution!\n');
        break;
    end
    w=w-alpha*(dfw+v+beta/2*(w-W*w));

    v=v+beta/2*(w-W*w);

    if strcmp(loss,'l2')
        Error_Q(i_iteration)=norm(w-threshold)^2/n;
    elseif strcmp(loss,'l1')
        Error_Q(i_iteration)=norm(w-threshold,1)/n;
    elseif strcmp(loss,'inf')
        Error_Q(i_iteration)=max(abs(w-threshold));
    else
        fprintf('Please input right loss type! Options include l2, l1 and inf! \n');
        break;
    end
    % max_value=max(abs(w-threshold));

    i_iteration=i_iteration+1
end