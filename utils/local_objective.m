function rho=local_objective(x,p)

L=length(x);

rho=zeros(L,1);

for i=1:L
    if x(i)<=0
        rho(i)=(p-1)*x(i);
    else
        rho(i)=p*x(i);
    end
end