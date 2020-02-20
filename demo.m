clear all

if 1 %Gaussian dataset
    p=300; n=5000;                         % p, number of samples; n, number of variables
    rng('default');
    D=randn(p,n);                       % data matrix
    D=D-repmat((mean(D,1)),p,1);        % Centering of the data
    A=D'*D;
    [d,ix]=sort(diag(A),'descend');A=A(ix,ix);D=D(:,ix);
end

result=[];
%set parameters for TPower
options.cardinality_vec = 100;
options.OptTol = 1e-6;
options.MaxIter = int32(500);
% run our method
tic
[U,f_result] = our_method(A, options,int32(500),true);
toc
x=U(:,end); %solution
hold on
plot(1:length(f_result),f_result,'m');
ylabel('Objective');
xlabel('Iteration');
plot(1:length(f_result),linspace(f_result(1),f_result(1),length(f_result)),'k--');
legend('Ours','TPower','Location','East');
title('Objective vs Iteration');