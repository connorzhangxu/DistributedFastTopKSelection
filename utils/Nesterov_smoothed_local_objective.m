function srho=Nesterov_smoothed_local_objective(x,p,mu)

L=length(x);

srho=zeros(L,1);

for i=1:L
  
    if x(i) <mu*(p-1)
        srho(i)=(p-1)*x(i)-mu/2*(p-1)^2;
    elseif x(i) > mu*p
        srho(i)=p*x(i)-mu/2*p^2;
    else
        srho(i)=x(i)^2/2/mu;
    end
    
end