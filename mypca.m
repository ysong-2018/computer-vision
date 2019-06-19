function [components variances] = pca(X,m)

% YOUR CODE HERE

%Center X, i.e., subtract the mean data point from each data point (in other words,
%subtract from each row the mean row in your data matrix). Check MATLAB's handy
%mean function.

meanrow = mean(X);
N = size(X,1);
X_C = zeros(size(X,1),size(X,2));

for i=1:N
    X_C(i,:) = X(i,:) - meanrow;
end

%Construct the covariance matrix C = 1/(n-1) XTX.
C = transpose(X_C)*X_C  / (N-1);

%Compute its eigenvectors and eigenvalues, i.e., the principal components and principal
%values of the original matrix X. MATLAB's eigs function allows you to specify how
%many eigenvectors/eigenvalues you want.
%returns diagonal matrix D of eigenvalues and matrix V whose columns are the corresponding right eigenvectors
[V,D] = eigs(C,m);
components = V;


%Convert the matrix of principal values returned in 3 into a vector. To do this, ex-
%tract the diagonal with diag. These are the variances of each corresponding principal
%component.
variances = diag(D);

end