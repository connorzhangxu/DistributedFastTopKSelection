%% Assign N numbers to n nodes, each node has more than 1 number
function  data=RandAssignData(N,n,x)
    random_numbers = randperm(N-1, n-1);
    sorted_numbers = sort(random_numbers);
    data={};
    data{1}=x(1:sorted_numbers(1));
    for i=2:n-1
        data{i}=x(sorted_numbers(i-1)+1:sorted_numbers(i));
    end
    data{n}=x(sorted_numbers(n-1)+1:N);
end