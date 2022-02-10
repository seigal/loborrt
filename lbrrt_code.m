% matlab code for
% "Lower bounds on the rank and symmetric rank of real tensors"
% by Kexin Wang and Anna Seigal

% to illustrate the results in Section 5,
% to find a real counter-example to Comon's conjecture of order 6

% we use the tensorlab toolbox for their functions cpdgen and tens2mat

%% Figure 2: illustrating the sets W1 and W2

% set the value of n
n = 6;

% construct alpha_i in R^{2n} as the columns of matrix alpha
alpha = zeros(2*n,n); 
for i = 1:n
    alpha(2*i-1,i) = 1;
    alpha(2*i,i) = 1;
end

% initialising the 16 types of vector u with u^{\otimes 5} in W1 
W1 = {};
for i = 1:16
    W1{i} = zeros(4*n,0);
end

% constructing the vectors u with u^{\otimes 5} in W1
% vectors u depending on four i-indices
t1 = 1;
t5 = 1;
for i1 = 1:n
    for i2 = i1+1:n 
        for i3 =i2+1:n
            for i4 = i3+1:n
                b = zeros(4*n,1);
                b(1:2*n) = alpha(:,i1) + alpha(:,i2) + alpha(:,i3) + alpha(:,i4);
                W1{5}(:,t5) = b;
                t5 = t5+ 1;
                for k1 = 1:n
                    b = zeros(4*n,1);
                    b(1:2*n) = alpha(:,i1) + alpha(:,i2) + alpha(:,i3) + alpha(:,i4);
                    b(2*n+1:4*n) = alpha(:,k1);
                    W1{1}(:,t1) = b;
                    t1 = t1+ 1;
                end
            end
        end
    end
end
% vectors u depending on three i-indices
t2 = 1;
t6 = 1;
t10 = 1;
t13 = 1;
for i1 = 1:n
    for i2 = i1+1:n
        for i3 =i2+1:n
            b = zeros(4*n,1);
            b(1:2*n) = alpha(:,i1) + alpha(:,i2) + alpha(:,i3);
            W1{6}(:,t6) = b;
            t6 = t6 + 1;
            for k1 = 1:n
                b = zeros(4*n,1);
                b(1:2*n) = alpha(:,i1) + alpha(:,i2) + alpha(:,i3);
                b(2*n+1:4*n) = ((n-3)/(n-4))*alpha(:,k1);
                W1{2}(:,t2) = b;
                t2 = t2+ 1;
                b(2*n+1:4*n) = ((n-2)/(n-1))*alpha(:,k1);
                W1{13}(:,t13) = b;
                t13 = t13 + 1;
                for k2 = k1+1:n
                    b = zeros(4*n,1);
                    b(1:2*n) = alpha(:,i1) + alpha(:,i2) + alpha(:,i3);
                    b(2*n+1:4*n) = alpha(:,k1) + alpha(:,k2);
                    W1{10}(:,t10) = b;
                    t10 = t10 + 1;
                end
            end
        end
    end
end
% vectors u depending on two i-indices
t11 = 1;
t3 = 1;
t7 = 1;
t14 = 1;
for i1 = 1:n
    for i2 = i1+1:n
        b = zeros(4*n,1);
        b(1:2*n) = alpha(:,i1) + alpha(:,i2);
        W1{7}(:,t7) = b;
        t7 = t7 + 1;
        for k1 = 1:n
            b = zeros(4*n,1);
            b(1:2*n) = alpha(:,i1) + alpha(:,i2);
            b(2*n+1:4*n) = ((n-2)/(n-4))*alpha(:,k1);
            W1{3}(:,t3) = b;
            t3 = t3 + 1;
            b(2*n+1:4*n) = ((n-2)^2)/((n-3)*(n-1))*alpha(:,k1);
            W1{14}(:,t14) = b;
            t14 = t14 + 1;
            for k2 = k1+1:n
                b = zeros(4*n,1);
                b(1:2*n) = alpha(:,i1) + alpha(:,i2);
                b(2*n+1:4*n) = ((n-2)/(n-3))*(alpha(:,k1) + alpha(:,k2));
                W1{11}(:,t11) = b;
                t11 = t11 + 1;
            end
        end
    end
end
% vectors u depending on one i-index
t4 = 1;
t8 = 1;
t12 = 1;
t15 = 1;
for i1 = 1:n
    b = zeros(4*n,1);
    b(1:2*n) = alpha(:,i1);
    W1{8}(:,t8) = b;
    t8 = t8 + 1;
    for k1 = 1:n
        b = zeros(4*n,1);
        b(1:2*n) = alpha(:,i1);
        b(2*n+1:4*n) = ((n-1)/(n-4))*alpha(:,k1);
        W1{4}(:,t4) = b;
        t4 = t4 + 1;
        b(2*n+1:4*n) = ((n-2)/(n-3))*alpha(:,k1);
        W1{15}(:,t15) = b;
        t15 = t15 + 1;
        for k2 = k1+1:n
            b = zeros(4*n,1);
            b(1:2*n) = alpha(:,i1);
            b(2*n+1:4*n) = ((n-1)/(n-3))*(alpha(:,k1) + alpha(:,k2));
            W1{12}(:,t12) = b;
            t12 = t12 + 1;
        end
    end
end
% vectors u depending on no i-indices
t9 = 1;
t16 = 1;
for k1 = 1:n
    b = zeros(4*n,1);
    b(2*n+1:4*n) = alpha(:,k1);
    W1{9}(:,t9) = b;
    t9 = t9 + 1;
    for k2 = k1+1:n
        b = zeros(4*n,1);
        b(2*n+1:4*n) = alpha(:,k1) + alpha(:,k2);
        W1{16}(:,t16) = b;
        t16 = t16 + 1;
    end
end

% collecting the u vectors together into the matrix W1all 
W1all = zeros(4*n,0);
for i = 1:16
    W1all = [W1all, W1{i}];
end

% constructing the vectors u with u^{\otimes 5} in W2 from the matrix W1all 
W2all = zeros(size(W1all));
W2all(1:2*n-1,:) = W1all(2*n+2:end,:);
W2all(2*n,:) = W1all(2*n+1,:);
W2all(2*n+1:4*n-1,:) = W1all(2:2*n,:);
W2all(4*n,:) = W1all(1,:);

% plotting the matrices W1all and W2all
subplot(2,1,1)
imagesc(W1all)
axis off
subplot(2,1,2)
imagesc(W2all)
axis off

%% Proposition 5.3: the clones of certain monomials are in span W

%% Lemma A.1
% the clone of x^4 y is in span W1
% note: uses tensorlab for the function cpdgen
% note: uses the list W1, see above, for some n \geq 5

% building 16 tensors sum u^{\otimes 5}, for the 16 types of W1 vector
T = {};
for i = 1:16
    T{i} = cpdgen({W1{i},W1{i},W1{i},W1{i},W1{i}});
end

% summing them together to get the linear combination 
coeff{1} = 1;
coeff{2} = -((n-4)^2)/(n-3);
coeff{3} = ((n-3)*(n-4)^2)/(2*(n-2));
coeff{4} = -((n-2)*(n-3)*(n-4)^2)/(6*(n-1));
coeff{5} = -n;
coeff{6} = (n-4)^2*n/(n-3);
coeff{7} = -(n-3)*(n-4)^2*n/(2*(n-2));
coeff{8} = ((n-2)*(n-3)*((n-4)^2)*n)/(6*(n-1));
coeff{9} = -(nchoosek(n,4) - ((n-3)^4/(n-4)^3)*nchoosek(n,3) + ((n-3)*(n-2)^4/(2*(n-4)^3))*nchoosek(n,2) - (n-2)*(n-3)*((n-1)^4)*n/(6*(n-4)^3));
bigT = zeros(4*n,4*n,4*n,4*n,4*n);
for i = 1:9
    bigT = bigT + coeff{i}*T{i};
end

% showing that the tensor bigT is a clone, by constructing 2x2x2x2x2 tensors
% of the minimum and maximum on each 2nx2nx2nx2nx2n block, and showing that 
% the two tensors are the same
E = [1:2*n; 2*n+1:4*n]; % first row is indices in E second row is indices in mathcal(E)
minT = zeros(2,2,2,2,2);
maxT = zeros(2,2,2,2,2);
for i1 = 1:2
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                for i5 = 1:2
                    subT = bigT(E(i1),E(i2),E(i3),E(i4),E(i5));
                    subT = subT(:);
                    minT(i1,i2,i3,i4,i5) = min(subT);
                    maxT(i1,i2,i3,i4,i5) = max(subT);
                end
            end
        end
    end
end
norm(minT(:)-maxT(:)) % this norm is zero, so bigT is a clone 
% plotting the flattening of minT, to see that it is x^4 y
minM = tens2mat(minT,[1 2],[3 4 5]);
imagesc(minM)

%% Lemma A.2

% the clone of x^3 y^2 is in span W1
% note: uses tensorlab for the function cpdgen
% note: uses the list W1, see above, for some n \geq 4

% building 16 tensors sum u^{\otimes 5}, for the 16 types of W1 vector
S = {};
for i = 1:16
    S{i} = cpdgen({W1{i},W1{i},W1{i},W1{i},W1{i}});
end
coeff{1} = 0;
coeff{2} = 0;
coeff{3} = 0;
coeff{4} = 0;
coeff{5} = 0;
coeff{6} = -(nchoosek(n,2)-(n-1)^2*n/(n-2));
coeff{7} = -(-nchoosek(n,2)*(n-3)^3/(n-2)^2+(n-1)^2*(n-3)^3*n/(n-2)^3);
coeff{8} = -(nchoosek(n,2)*(n-2)*(n-3)^3/(2*(n-1)^2)-(n-3)^3*n/2);
coeff{9} = -(-((n-2)^4/(n-1)^3)*nchoosek(n,3)+((n-2)^7/((n-3)^2*(n-1)^3))*nchoosek(n,2)-(n-2)^5*n/(2*(n-3)^2));
coeff{10} = 1;
coeff{11} = -(n-3)^3/(n-2)^2;
coeff{12} = (n-2)*(n-3)^3/(2*(n-1)^2);
coeff{13} = -(n-1)^2/(n-2);
coeff{14} = (n-1)^2*(n-3)^3/(n-2)^3;
coeff{15} = -(n-3)^3/2;
coeff{16} = -(nchoosek(n,3)-nchoosek(n,2)*((n-2)^3/(n-3)^2)+(n-1)^3*(n-2)*n/(2*(n-3)^2));

% summing them together to get the linear combination 
bigS = zeros(4*n,4*n,4*n,4*n,4*n);
for i = 1:16
    bigS = bigS + coeff{i}*S{i};
end

% showing that the tensor bigS is a clone 
E = [1:2*n; 2*n+1:4*n];
minS = zeros(2,2,2,2,2);
maxS = zeros(2,2,2,2,2);
for i1 = 1:2
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                for i5 = 1:2
                    subS = bigS(E(i1),E(i2),E(i3),E(i4),E(i5));
                    subS = subS(:);
                    minS(i1,i2,i3,i4,i5) = min(subS);
                    maxS(i1,i2,i3,i4,i5) = max(subS);
                end
            end
        end
    end
end
norm(minS(:)-maxS(:)) % this is zero, so bigS is a clone 
minN = tens2mat(minS,[1 2],[3 4 5]);
imagesc(minN)


%% Figure 3: the set W^{(4)} is linearly independent

% construct the vectors from Shitov's equation (11.1)
n = 5;
beta=zeros(10,n);
for i = 1:5
    beta(2*i-1,i) = 1;
    beta(2*i,i) = 1;
end

% initialising the 5 types of vectors u with u^{\otimes 3} in W1 
W1 = {};
for i = 1:16
    W1{i} = zeros(4*n,0);
end

% vectors depending on three E-indices
t1 = 1;
t2 = 1;
for i = 1:n
    for j = i+1:n
        b = zeros(20,1);
        b(1:10) = beta(:,i) + beta(:,j);
        W1{2}(:,t2) = b;
        t2 = t2 + 1;
        for k = 1:5
            b = zeros(20,1);
            b(1:10) = beta(:,i) + beta(:,j);
            b(11:20) = beta(:,k);
            W1{1}(:,t1) = b;
            t1 = t1 + 1;
        end
    end
end
% vectors depending on three E-indices
t3 = 1;
t4 = 1;
for i = 1:n
    b = zeros(20,1);
    b(1:10) = beta(:,i);
    W1{4}(:,t4) = b;
    t4 = t4 + 1;
    for k = 1:n
        b = zeros(20,1);
        b(1:10) = 3*beta(:,i);
        b(11:20) = 4*beta(:,k);
        W1{3}(:,t3) = b;
        t3 = t3 + 1;
    end
end
% vectors depending on one E-index
t5 = 1;
for i = 1:n
    b = zeros(20,1);
    b(11:20) = beta(:,i);
    W1{5}(:,t5) = b;
    t5 = t5 + 1;
end

W1all = zeros(4*n,0);
for i = 1:5;
    W1all = [W1all, W1{i}];
end

% the vectors u^{\otimes 3}, for u \in span W1, are the columns of the 
% 8000 x 95 matrix M 
ncol = size(W1all,2);
M = zeros(power(20,3),ncol);
for i = 1:ncol
    v = cpdgen({W1all(:,i),W1all(:,i),W1all(:,i)});
    M(:,i) = v(:);
end

Moriginal = M;
Msupp = (Moriginal > 0);
Msum={};
Msum{1} = sum(Msupp,2);
i=1;

while 1
    % removing the entries that are one coefficient 
    k = find(Msum{i} == 1);
    length(k); % the number of entries that are one coefficient 
    if isempty(k) % everything is zero so we get linear independence 
        i; % number of loops
        break
    end
    Mnew = Moriginal;
    colinds = zeros(length(k),1);
    for j = 1:length(k)
       Mk = Msupp(k(j),:);
       colind = find(Mk > 0);
       colinds(j) = colind;
       Mnew(:,colind) = zeros(8000,1);
    end
    Moriginal = Mnew;
    length(unique(colinds)); % number of rank one tensors set to 0 in this step
    clear Msupp
    i = i+1;
    Msupp = (Mnew > 0);
    Msum{i} = sum(Msupp,2);
end


% plotting the supports of the matrices in the 3 steps as 20x400 matrices
% plotting log(1 + B) to see support clearly 
C1 = log(reshape(Msum{1},20,400) + ones(20,400));
C2 = log(reshape(Msum{2},20,400) + ones(20,400));
C3 = log(reshape(Msum{3},20,400) + ones(20,400));
C4 = log(reshape(Msum{4},20,400) + ones(20,400));
clims = [0 2];
subplot(3,1,1)
imagesc(C1,clims)
axis off
subplot(3,1,2)
imagesc(C2,clims)
axis off
subplot(3,1,3)
imagesc(C3,clims)
axis off




