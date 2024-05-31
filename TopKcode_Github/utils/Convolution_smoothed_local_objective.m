function lh=Convolution_smoothed_local_objective(x,p,h)

L=length(x);

lh=zeros(L,1);

for i=1:L
  
    if x(i) <= -h 
        lh(i)=(p-1)*x(i);
    elseif x(i) <= h
        lh(i)= (x(i)-h)^2/4/h+p*x(i);
    else
        lh(i)=p*x(i);
    end
    
end