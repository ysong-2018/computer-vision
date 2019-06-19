%% Homework 4
% * Name: Audrey Song
% * NetId: ys585

%IMPORTANT: Please do not modify any of the question headings!
%Also, please do not modify any of the file names (including this one).

%% 1a
% Your answer here
% I see that there is a distinction between the patterns. The first 9 are
% fairly similar, while the last 11 have a very distinct y value represented 
% by the yellow color. 

%Your code here
load tinypics;
imagesc(data);
axis equal tight;



%% 1b
% Your answer here
% I would like to say that there are 2 clusters. It is obvious when the 
% 3d plot is dragged.

%Your code here

x = data(:,1);
y = data(:,2);
z = data(:,3);

scatter3(x,y,z);
xlabel('pixel1');
ylabel('pixel2');
zlabel('pixel3');


%% 1c 
%Write your function using the mykmeans.m file. Please do not rename that file.



%% 1d 
% Your answer here
% I ran the data on mykmeans function and the results are plotted. The
% results basically match my expectations. Running it several times could
% give some different results, since in the algorithm, clustering happens
% by starting at a random point, that could potentially give different
% answer. 

%Your code here
[clusters, centroids] = mykmeans(data, 2);

firstcluster = data(clusters == 1,:);
x1 = firstcluster(:,1);
y1 = firstcluster(:,2);
z1 = firstcluster(:,3);
scatter3(x1,y1,z1,20,'g');
hold on;

secondcluster = data(clusters == 2,:);
x2 = secondcluster(:,1);
y2 = secondcluster(:,2);
z2 = secondcluster(:,3);

scatter3(x2,y2,z2,20,'b');
hold on;

scatter3(centroids(:,1),centroids(:,2),centroids(:,3), 30,'r');

%%
close all; clear all

%% 2a (optional)
% Your answer here
% Here that 10 images are plotted. Using my brain, I would like to separate
% it into two classes. The first five are more clear and the last five are
% kind of blurry and have a bigger Y.

%Your code here
load logos;

figure
subplot(2,5,1);imshow(data(:,:,1));
subplot(2,5,2);imshow(data(:,:,2));
subplot(2,5,3);imshow(data(:,:,3));
subplot(2,5,4);imshow(data(:,:,4));
subplot(2,5,5);imshow(data(:,:,5));
subplot(2,5,6);imshow(data(:,:,6));
subplot(2,5,7);imshow(data(:,:,7));
subplot(2,5,8);imshow(data(:,:,8));
subplot(2,5,9);imshow(data(:,:,9));
subplot(2,5,10);imshow(data(:,:,10));


%% 2b (optional)
% Your answer here
% Here the X shown is the data matrix required. The dimensionality for each
% row in my matrix is 1 by 4800.
% I am not able to make the data into a scatterplot like 1-b, because each
% image has too many data values and i cannot plot a 4800 dimensional
% scatterplot

%Your code here
X = zeros(10,4800);

for i=1:10
    cur = data(:,:,i);
    vec_cur = transpose(cur(:));
    X(i,:) = vec_cur;
end


%% 2c (optional)
% Your answer here
% I followed the instructions by running kmeans on the matrix. The clusters
% results given in idx sometimes don't match my expectations but after the
% kmeans has been run a few times, the results match my expectations
% better.
% The two centroids are plotted. 

%Your code here

[idx,C] = kmeans(X,2,'Start','uniform');
centroid1 = reshape(C(1,:),[48,100]);
centroid2 = reshape(C(2,:),[48,100]);


figure
subplot(121); imshow(centroid1);
subplot(122); imshow(centroid2);

%%
close all; clear all

%% 3a
%Write your function using the the mypca.m file. Please do not rename that file.

%% 3b
% Your answer here
% All the required are implemented

%Your code here
load gaussian;

[components, variances] = mypca(gaussian , 2);

scaled_components = components;
scaled_components(:,1) = components(:,1)*variances(1);
scaled_components(:,2) = components(:,2)*variances(2);
x=[0 0];
y=[0 0];

scatter(gaussian(:,1),gaussian(:,2));
hold on;
quiver(x,y,scaled_components(1,:),scaled_components(2,:));
axis equal;

%% 3c
% Your answer here
% 
% Covariance provides a measure of the strength of the correlation between
% the vectors. 
% A covariance matrix is a matrix whose element in the i, j position is 
% the covariance between the i-th and j-th elements of a random vector.
% it is important that each component of the vectors in the data matrix 
% have zero mean because we do not want the data to be all too large or all
% too small. It is best to have it centered at zero. 

%% 3d
% Your answer here
% The related singular value decomposition values are calculated, and we
% can see that the components and variances are similar. This is because
% svd and pca are both dimension reduction methods. PCA uses svd in the
% calculations, but requires extra steps. 

%Your code here
[U,S,V] = svd(gaussian);
s = diag(S);
v = (s.^2)/(1000-1);


%%
close all; clear all

%% 4a
% Your answer here
% mypca has been run and and the four components are visualized. 
% I see that only the first component returned is quite significant.
% Also, only the first variance is large and other get small really
% quickly.
% I think a one dimensional array would suffice to represent the dataset,
% and the result is what I expected. 

%Your code here -- please reload the logos dataset for this part!
%(i.e. use a fresh variable for your data matrix, for grading purposes)

load logos;
logo = data;
X = zeros(10,4800);

for i=1:10
    cur = logo(:,:,i);
    vec_cur = transpose(cur(:));
    X(i,:) = vec_cur;
end


[components, variances] = mypca(X,4);
component1 = components(:,1);
component2 = components(:,2);
component3 = components(:,3);
component4 = components(:,4);
r1 = reshape(component1, [48,100]);
r2 = reshape(component2, [48,100]);
r3 = reshape(component3, [48,100]);
r4 = reshape(component4, [48,100]);


figure
subplot(221); imagesc(r1);
subplot(222); imagesc(r2);
subplot(223); imagesc(r3);
subplot(224); imagesc(r4);

figure
plot(variances);


%% 4b
% Your answer here
% After projecting, the dimensions of the data becomes 10 by 1. 
% The relation between the components matrix and the concept of basis or
% coordinate system: we can view the columns of the components as the
% coordinate system that we can manipulate on. 

%Your code here
projection = X*component1;



%% 4c
% Your answer here
% The scatter plot is done and we can see two clusters. The resulting 
% assignment what I expected based on what the images look like. 
% The average images are also plotted. After running a few times, we always
% see that the first centroid is representative of one cluster and
% the second centroid is representative of the other. 

%Your code here
[idx,C] = kmeans(projection,2,'Start','uniform');


%scatter(idx,zrs);

firstcluster = projection(idx == 1);
zrs1 = zeros(size(firstcluster,1),1);
scatter(firstcluster,zrs1,20,'g');
hold on;

secondcluster = projection(idx == 2);
zrs2 = zeros(size(secondcluster,1),1);
scatter(secondcluster,zrs2,20,'b');
hold on;

scatter(C,[0 0], 35,'r');


first = X(idx == 1,:);
second = X(idx == 2,:);

first = mean(first);
second = mean(second);

i1 = reshape(first, [48,100]);
i2 = reshape(second, [48,100]);

figure
subplot(121);imshow(i1);
subplot(122);imshow(i2);



%% 4d
% Your answer here
% The dimensionality reduction works better. This is because in higher
% dimensions, distance measure does not work very well, making clustering
% less accurate. However, PCA does not have this issue and things can be
% done easily through matrix manipulations. 

%%
close all; clear all

%% 5a
% Your answer here
% Instructions are followed and eigenfaces are plotted.
% The eigenfaces tell me the most significant features of people's faces,
% and the relevant features when comparing face images are the contours of
% faces, locations of eyes, nose, mouths, and other facial features. 

%Your code here
load faces;
face = double(Data);

for i =1:size(face,1)
    norm = max(face(i,:));
    disp(norm);
    face(i,:) = face(i,:)/norm;
end

imshow(transpose(reshape(face(23,:), [96,96])));


[components, variances] = mypca(face,9);

component1 = components(:,1);
component2 = components(:,2);
component3 = components(:,3);
component4 = components(:,4);
component5 = components(:,5);
component6 = components(:,6);
component7 = components(:,7);
component8 = components(:,8);
component9 = components(:,9);
r1 = transpose(reshape(component1, [96,96]));
r2 = transpose(reshape(component2, [96,96]));
r3 = transpose(reshape(component3, [96,96]));
r4 = transpose(reshape(component4, [96,96]));
r5 = transpose(reshape(component5, [96,96]));
r6 = transpose(reshape(component6, [96,96]));
r7 = transpose(reshape(component7, [96,96]));
r8 = transpose(reshape(component8, [96,96]));
r9 = transpose(reshape(component9, [96,96]));

figure
subplot(3,3,1); imshow(r1,[]);
subplot(3,3,2); imshow(r2,[]);
subplot(3,3,3); imshow(r3,[]);
subplot(3,3,4); imshow(r4,[]);
subplot(3,3,5); imshow(r5,[]);
subplot(3,3,6); imshow(r6,[]);
subplot(3,3,7); imshow(r7,[]);
subplot(3,3,8); imshow(r8,[]);
subplot(3,3,9); imshow(r9,[]);

figure
plot(variances);



%% 5b (optional)
% Your answer here
%Please write your function in a file named pcasvd.m. Please do not rename that file.

%Your code here

%%
close all; clear all


