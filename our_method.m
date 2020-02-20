function [X,f_result] = our_method(A, options,out_iter_num,warm_start)

%%   Enhancing Truncated Power Method for Sparse PCA by Support Altering
%     max x'*A*x    subject to ||x||=1, ||x||_0 <= k
%
% *** inputs ***
% - A:                      p x p covariance matrix
% - options:                a structure stores user specified parameters which include:
%    -- cardinality_vec:    an m-dimensional vector stores the cadinality for each sparse loading
%    -- optTol:             optimality tolerance
%    -- maxIter:            max number of iteration in each TPower
% - out_iter_num:           iteration limit of our algorithm
% - warm_start:             warm start strategy (0: do not employ warm start; 1: employ warm start)
% 
% *** outputs ***
% - X:                      each column of X is a sparse vector
% - f_result:               each column of f_result is the objective value corresponding to the sparse vector

%% Set Parameters
cardinality_vec=options.cardinality_vec;
optTol=options.OptTol;
maxIter=options.MaxIter;

%% Initialization
X = [];
dim = size(A,1);
%% Main loop to extract multiple sparse loadings
cardinality = cardinality_vec;
options_cur.cardinality = cardinality;
options_cur.optTol = optTol;
options_cur.maxIter = maxIter;

%warm start selection
if ~warm_start
    [x, f] = TPower_modified(A, options_cur);
else
    %with 8k 4k 2k k warm-start
    b=floor(dim/cardinality);
    if b>=8
        c=8;
        iteration=4;
    elseif b>=4
        c=4;
        iteration=3;
    elseif b>=2
        c=2;
        iteration=2;
    else
        c=1;
        iteration=1;
    end
    if c~=8
        [x,f]=eigs(A,1);
    else
        options_cur.cardinality = c*cardinality;
        [x, f] = TPower_modified(A, options_cur);
    end
    for i=1:iteration
        options_cur.cardinality = c*cardinality;
        [x, f] = TPower_modified(A, options_cur,x);
        c=c/2;
    end
end

x_init=x;
a=min(cardinality,dim-cardinality);
change_num=a;
X = [X,x];
f_result=f;
max_changes=a;
for i=1:out_iter_num-1
    j=1;
    [x2,out_X]=sa_partial(x_init,A,max_changes);
    while j<=max_changes
        [x, f] = TPower_modified(A, options_cur,x2);
        %Is the value higher?
        if f>f_result(i)
            %save result
            f_result=[f_result;f];
            x_init=x;
            X = [X,x];
            max_changes=change_num;%non-strictly decreasing
            break;
        else
            %decrease the change_num
            change_num=max(max_changes-j,1);
            x2=out_X(:,change_num);
            j=j+1;
        end
    end
    if j>max_changes
        break;
    end
end
%make sure cw-max
while true
    x2=sa_cwmax(x_init,A);
    [x, f] = TPower_modified(A, options_cur,x2);
    %Is the value higher?
    if f>f_result(i)
        %save result
        f_result=[f_result;f];
        x_init=x;
        X = [X,x];
        i=i+1;
    else
        break;
    end
end
x=x_init;%final sparse pc loading
end
