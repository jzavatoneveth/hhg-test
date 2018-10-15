function [p, t, pstat] = HHGPermutationTest(X,Y,nperm,maxN)
%HHGPermutationTest Heller-Heller-Gorfine multivariate test of association
%
%   Syntaxes:
%       [ p ] = HHGPermutationTest(X,Y) computes the p-value of the HHG
%           test using 100 permutations for an (Nxm1) matrix X and an
%           (Nxm2) matrix Y.
%
%       [ p ] = HHGPermutationTest(X,Y,nperm) computes the p-value of the
%           HHG test using nperm permutations
%
%       [ p ] = HHGPermutationTest(X,Y,nperm,maxN) computes the p-value of
%           the HHG test while storing at most a (maxN x maxN) distance
%           matrix in memory. If n > maxN, incremental computation is used.
%
%       [p, t] = HHGPermutationTest(...);
%       [p, t, pstat] = HHGPermutationTest(...);
%
%   Inputs:
%       'X', 'Y' - Input data matrices. Rows are assumed to represent
%           samples, and columns are assumed to represent dimensions.
%
%       'nperm' - The number of permutations to use in computing the exact
%           test. The default value is 100.
%
%       'maxN' - The maximum number of input samples for which the entire
%           distance matrix will be stored in memory. The default value is
%           1e4.
%
%   Outputs:
%       'p' - The p-value of the permutation test.
%
%       't' - The HHG test statistic.
%
%       'pstat' - A nperm by 1 vector of the HHG test statistics for each
%           permutation.
%
%   References:
%       [1] - Heller, R., Heller, Y., & Gorfine, M. (2012). A consistent
%           multivariate test of association based on ranks of distances.
%           Biometrika, 100(2), 503-510.
%
%   Copyright (c) 2018 Jacob Zavatone-Veth, MIT License

% Set the  number of permutations
if (nargin < 3) || isempty(nperm)
    nperm = 100;
elseif (~isscalar(nperm)) || (rem(nperm,1) > 0) || (nperm < 1)
    error('The number of permutations must be a scalar integer greater than or equal to 1.');
end

% Set the maximum sample size for which the full NxN distance matrix will
% be held in memory 
if (nargin < 4) || isempty(maxN)
    maxN = 1e4;
elseif (~isscalar(maxN)) || (rem(maxN,1) > 0) || (maxN < 1)
    error('The maximum number of points for in-memory computation must be a scalar integer greater than or equal to 1.');
end

% Check that the input matrices contain the same number of samples
if (size(X,1) ~= size(Y,1))
    error('Input matrices must contain the same number of samples.');
end

% Get the number of samples
n = size(X,1);

% Get the field width of the number of permutations
ndig = 1+floor(log10(nperm));

% Allocate a container to store the values of the test statistic for each permutation
pstat = nan(nperm,1);

% Start a timer
tAll = tic;

% Decide whether to store the full pairwise distance matrices in memory
if n < maxN
    
    % Compute the full distance matrices
    fprintf('%d points: computing full pairwise distance matrices.\n', n);
    dx = squareform(pdist(X));
    dy = squareform(pdist(Y));
    fprintf('Computed distance matrices: %f seconds elapsed.\n', toc(tAll));
    
    % Compute the test statistic
    t = hhgTestStatisticFromFullDistanceMatrices(dx,dy);
    fprintf('Computed test statistic: %f seconds elapsed.\n', toc(tAll));
    
    % Compute the distribution of the test statistic under permutation of y
    for ind = 1:nperm
        idx = randperm(n);
        pstat(ind) = hhgTestStatisticFromFullDistanceMatrices(dx, dy(idx,idx));
        fprintf('Permutation %*d of %d: %f seconds elapsed.\n', ndig, ind, nperm, toc(tAll));
    end
    
else
    
    % If we fall into this case, the distance matrix is too large, and we
    % must compute the test statistic incrementally.
    fprintf('%d points: using incremental algorithm.\n', n);
    
    % Incrementally compute distance matrices and write them to files
    [ posx, posy ] = hhgDistanceMatrixToFile(X,Y);
    fprintf('Computed distance matrices: %f seconds elapsed.\n', toc(tAll));
    
    % Compute the test statistic
    t = hhgTestStatisticIncremental(posx,posy,(1:n)');
    fprintf('Computed test statistic: %f seconds elapsed.\n', toc(tAll));
    
    % Compute the distribution of the test statistic under permutation of y
    for ind = 1:nperm
        idx = randperm(n);
        pstat(ind) = hhgTestStatisticIncremental(posx,posy,idx);
        fprintf('Permutation %*d of %d: %f seconds elapsed.\n', ndig, ind, nperm, toc(tAll));
    end
    
    % Delete temporary files
    delete('dx.dat');
    delete('dy.dat');
    
end

% Compute the p-value
p = nnz(pstat >= t) / nperm;

% Print a timing message to the terminal
fprintf('Computation finished in %f seconds\n', toc(tAll));

end

function [ t ] = hhgTestStatisticFromFullDistanceMatrices(dx,dy)
% Local utility function to compute the HHG test statistic from precomputed
% distance matrices

% Get the number of samples
n = size(dx,1);
m = n - 1;

% Compute needed indexing vectors
c = (1:n)';
j = (1:m)';

% Initialize the test statistic
t = 0;

% Iterate over samples
for i = 1:n
    % Extract the desired data from the distance matrices
    dx_tmp = dx(c~=i, i);
    dy_tmp = dy(c~=i, i);
    
    % Re-order the distances by the distance to i in x
    [dx_tmp, jIdx] = sort(dx_tmp, 'ascend');
    dy_tmp = dy_tmp(jIdx);
    
    % Compute the permutation pi(j) defined as the rank of the distances to i in y
    [~,pij] = ismember(dy_tmp, sort(dy_tmp, 'ascend'));
    
    % Count the number of inversions for each j
    invj = double(inversions(uint64(pij)));
    
    % Compute the contingency table
    a11 = j - 1 - invj;
    a12 = invj;
    a21 = pij + invj - j;
    a22 = n - 1 - pij - invj;
    
    % Compute the classical test statistic
    snum = (a12 .* a21 - a11 .* a22).^2;
    sden = (j-1).*(n-j-1).*(pij-1).*(n-1-pij);
    s = snum ./ sden;
    s(sden==0) = 0;
    
    % Accumulate the test statistic
    t = t + sum(s);
end

% Adjust the test statistic
t = (n-2)*t;

end

function [ posx, posy ] = hhgDistanceMatrixToFile(x,y)
% Incrementally compute distance matrices and write them to disk

% If distance matrix files already exist, delete them
if (exist('dx.dat', 'file') == 2), delete('dx.dat'); end
if (exist('dy.dat', 'file') == 2), delete('dy.dat'); end

% Open file IDs
fx = fopen('dx.dat','w');
fy = fopen('dy.dat','w');

% Get the number of samples
n = size(x,1);

% Allocate containers to store offsets
posx = nan(n,1);
posy = nan(n,1);

% Iterate over samples
for i = 1:n
    % Get current offsets
    posx(i) = ftell(fx);
    posy(i) = ftell(fy);
    
    % Compute distances in x and y from point i
    dx = pdist2(x, x(i,:));
    dy = pdist2(y, y(i,:));
    
    % Write distances to files
    fwrite(fx, dx, 'double');
    fwrite(fy, dy, 'double');
    
end

% Close file IDs
fclose(fx);
fclose(fy);

end


function [ t ] = hhgTestStatisticIncremental(posx, posy, idx)
% Local utility function to compute the HHG test statistic incrementally
% from distance matrices stored on disk

% Get the number of samples
n = size(posx,1);
m = n - 1;

% Compute needed indexing vectors
c = (1:n)';
j = (1:m)';

% Initialize the test statistic
t = 0;

% Open file IDs
fx = fopen('dx.dat','r');
fy = fopen('dy.dat','r');

% Permute posy
posy = posy(idx);

% Iterate over samples
for i = 1:n
    % Seek to the needed locations in the distance matrix files
    fseek(fx, posx(i), 'bof');
    fseek(fy, posy(i), 'bof');
    
    % Read in distance matrix data from file
    dx = fread(fx, n, 'double');
    dy = fread(fy, n, 'double');
    
    % Permute y-distances
    dy = dy(idx);
    
    % Remove the current point
    dx = dx(c~=i);
    dy = dy(c~=i);
    
    % Re-order the distances by the distance to i in x
    [dx, jIdx] = sort(dx, 'ascend');
    dy = dy(jIdx);
    
    % Compute the permutation pi(j) defined as the rank of the distances to i in y
    [~,pij] = ismember(dy, sort(dy, 'ascend'));
    
    % Count the number of inversions for each j
    invj = double(inversions(uint64(pij)));
    
    % Compute the contingency table
    a11 = j - 1 - invj;
    a12 = invj;
    a21 = pij + invj - j;
    a22 = n - 1 - pij - invj;
    
    % Compute the classical test statistic
    snum = (a12 .* a21 - a11 .* a22).^2;
    sden = (j-1).*(n-j-1).*(pij-1).*(n-1-pij);
    s = snum ./ sden;
    s(sden==0) = 0;
    
    % Accumulate the test statistic
    t = t + sum(s);
end

% Close file IDs
fclose(fx);
fclose(fy);

% Adjust the test statistic
t = (n-2)*t;

end
