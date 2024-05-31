% Generate random graph
function [A]=RandomGraphGeneration(N,NumEdges)
A=zeros(N);
diameter=N+1;
t=0;
while diameter>=N+1  % if the graph is not connected, the diameter will be Inf
    A=zeros(N);
    RealEdges=0;
    while RealEdges<NumEdges
        i=ceil(rand(1)*N);
        j=ceil(rand(1)*N);
        A(i,j)=1;
        A(j,i)=1;
        RealEdges=sum(sum(A))/2;
    end
    %% calculate the diameter
    DG= distances(graph(A));
    diameter=max(max(DG));
    t=t+1;
end
