function [clusters centroids] = mykmeans(X, k, init)
%MYKMEANS implements the k-means algorithm.
%   [CLUSTERS CENTROIDS] = MYKMEANS(X, K, INITS) partitions the data points
%   in the N-by-P data matrix X into K distinct clusters, using Euclidean
%   distance. This is a simple implementation of the k-means algorithm with
%   random initialization.
%
%   Optionally, it takes the argument INIT, a K-by-P matrix with a fixed
%   initial position for the cluster centroids.  
%
%   MYKMEANS returns an N-by-1 vector CLUSTERS containing the cluster
%   indices of each data point, as well as CENTROIDS, a K-by-P matrix with 
%   the final cluster centroids' locations.


[n p] = size(X);

if ~exist('init', 'var')
  %choose initial centroids by picking K points uniformly at random from the range of X
  init = min(X(:)) + rand(k,p)*(max(X(:))-min(X(:)));
end
%centroids is a k-by-p random matrix
%its i^th row contains the coordinates of the cluster with index i
centroids = init;

%initialize cluster assignment vector
clusters = zeros(n,1);

MAXITER = 1000;


for iter=1:MAXITER
    
    %create a new clusters vector to fill in with updated assignments
    new_clusters = zeros(n,1);
    

    %for each data point x_i
    for i=1:n
        
        x_i = X(i,:);
        
        %find closest cluster
        closest = findClosestCluster(x_i,centroids);%%%IMPLEMENT THIS FUNCTION AT THE END OF THIS FILE

        %reassign x_i to the index of the closest centroid found
        new_clusters(i) = closest;

    end
    
    if hasConverged(clusters,new_clusters)%%%IMPLEMENT THIS FUNCTION AT THE END OF THIS FILE
        %exit loop
        break 
    end
    
    %otherwise, update assignment
    clusters = new_clusters;
    %disp(clusters);
    %and recompute centroids
    centroids = recomputeCentroids(X,clusters,k);%%%IMPLEMENT THIS FUNCTION AT THE END OF THIS FILE

end

if iter == MAXITER
    disp('Maximum number of iterations reached!');
end

end

function closest = findClosestCluster(x_i,centroids)
% Compute Euclidean distance from x_i to each cluster centroid and return 
% the index of the closest one.

%%% Replace the following line with your own code
%closest = 1;

m = size(centroids,1);
smallest  = sqrt(sum((x_i-centroids(1,:)).^2));
idx = 1;
disp("disp(centroids(1,:));");
disp(centroids(1,:));

for j=2:m
    cur = sqrt(sum((x_i-centroids(j,:)).^2));
    if (cur < smallest)
        smallest = cur;
        idx=j;
    end
end

closest =idx;

end

function converged = hasConverged(old_assignment, new_assignment)
% Check if algorithm has converged, i.e., cluster assignments haven't
% changed since last iteration.

%%% Replace the following line with your own code

converged = true;
n = size(old_assignment,1);

for i = 1:n
    if old_assignment(i) ~= new_assignment(i)
        converged = false;
    end
end
    

end

function centroids = recomputeCentroids(X,clusters,k)
% Recompute centroids based on current cluster assignment.

%%% Replace the following line with your own code
 p=size(X,2);
 operations = zeros(k,p);
 
 numerator = zeros(k,p);
 denominator = zeros(k,1);
 
 N = size(clusters,1);
 disp(N);
 for i = 1:N
     group = clusters(i);
     numerator(group,:) = numerator(group,:) + X(i,:);
     denominator(group) = denominator(group) + 1;   
 end
 
 for j = 1:k
     operations(j,:) = numerator(j,:) / denominator(j);
 end

centroids = operations;


end
    