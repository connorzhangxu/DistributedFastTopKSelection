%% EXTRA
function [Error_Q]=DistributedQuantileEstimation_EXTRA(x,p,A,beta0,h,n_iteration,Delta,loss,smooth)

D=diag(sum(A));
L=D-A;
N=length(x);
[y,~]=sort(x,'ascend');
Error_Q=zeros(N,1);

v=zeros(N,1);
W=eye(N)-beta0*L;
w=x;
% beta=1/h/sqrt(1-Sigma(2));
beta=1/h;
alpha=1/beta;
max_value=inf;
i_iteration=1;

while i_iteration <=n_iteration %&& max_value> Delta/2

    if strcmp(smooth,'Nesterov')
        % Nesterov's smoothing
        dfw=((x-w)<=h*(p-1))*(1-p)-((x-w)>h*p)*p+(h*(p-1)<(x-w) & (x-w)<=h*p).*(w-x)/h;
    elseif strcmp(smooth,'Convolution')
        % Convolution smoothing
        dfw=((x-w)<=-h)*(1-p)-((x-w)>h)*p+(-h<(x-w) & (x-w)<=h).*((w-x+h)/2/h-p);
    else
        fprintf('Please input right smoothing type! Options include Nesterov and Convolution!\n');
        break;
    end

    w=w-alpha*(dfw+v+beta/2*(w-W*w));

    v=v+beta/2*(w-W*w);

    % w= (eye(N)+W)/2*w-alpha*dfw;
    %
    % dfw=((x-w)<=h*(p-1))*(1-p)-((x-w)>h*p)*p+(h*(p-1)<=(x-w) & (x-w)<=h*p).*(w-x)/h;
    %
    % w= W*x-alpha*dfw;

    if strcmp(loss,'l2')
        Error_Q(i_iteration)=norm(w-y(floor(p*N)+1))^2/N;
    elseif strcmp(loss,'l1')
        Error_Q(i_iteration)=norm(w-y(floor(p*N)+1),1)/N;
    elseif strcmp(loss,'inf')
        Error_Q(i_iteration)=max(abs(w-y(floor(p*N)+1)));
    else
        fprintf('Please input right loss type! Options include l2, l1 and inf! \n');
        break;
    end
    % max_value=max(abs(w-y(floor(p*N)+1)));

    i_iteration=i_iteration+1
end