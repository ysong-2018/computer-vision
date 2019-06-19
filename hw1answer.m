%% Homework 1
% * Name: Audrey Song
% * NetId: ys585

%IMPORTANT: Please do not modify any of the question headings!

%% 1a
% Your answer here:
% The two visualizations are shown below. Here are the differences: imshow
% displays how the original image looks like in grayscale;  imagesc
% visualizes by mapping the double precision numbers converted from im2double to a full range of colors in the colormap
% Note: Use "% " at the beginning of each line of text (don't forget the
% space).
%This is what happens if you don't use a space before your text.

%Your code here
im = imread('Paolina.tiff');
im = im2double(im);

figure
subplot(121); imshow(im);
subplot(122); imagesc(im);
axis image;



%% 1b
% Your answer here: 
% The contour is added to the image. What is displayed is a set of contour
% lines getting through the same "heights" of numbers in the matrix. In
% short, the contour function primarily displays some isolines of matrix
% im, where the isolines displayed are automatically chosen by matlab .


%Your code here
contour(im);axis ij;
axis image;


%% 1c
% Your answer here:
% The vector field is displayed using function quiver. After zooming on the
% image, I can see that there is a whole set of vectors sketched at each
% pixel. Then length of vectors represent the magnitude of the rate of
% change at a particular point, and the direction of the vector represent
% towards which direction the matrix numbers are changing. 

%Your code here
[Mx, My] = gradient(im);
quiver(Mx, My);axis ij;axis image;



%% 1d
% Your answer here: 
% The magnitute is calculated below. The magnitude represents the degree to
% which the matrix elements change, so when we visualized this, the parts
% where the matrix elements change more can be seen more apparently.

%Your code here

gmag = sqrt((Mx.^2)+(My.^2));
imagesc(gmag);axis image;



%%
close all



%% 2a
% Your answer here:
% The three visualizations for functions 1 and 2 are shown below. The main
% difference is that function 1 has a gradual change and one high hill,
% while function 2 is generally flat but has two peaks.

%Your code here
[HISTM, FOOD] =    acoli_hist([0,0], 1, 0.5, 0.5, 50, 300, 50);
figure
subplot(131); imagesc(FOOD);
subplot(132); contour(FOOD);axis ij;
subplot(133); surf(FOOD);axis ij;

[HISTM, FOOD] =    acoli_hist([0,0], 2, 0.5, 0.5, 50, 300, 50);
figure
subplot(131); imagesc(FOOD);
subplot(132); contour(FOOD);axis ij;
subplot(133); surf(FOOD);axis ij;

%% 2b
% Your answer here:
% The gradient vector fields are displayed. We can see that the gradient
% fields match the 3d visulizations, in that the gradient vectors have 
% greater magnitudes at the points where the food function experiences 
% faster change. 

%Your code here
[HISTM, FOOD] =    acoli_hist([0,0], 1, 0.5, 0.5, 50, 300, 50);
[Mx1, My1] = gradient(FOOD);

[HISTM, FOOD] =    acoli_hist([0,0], 2, 0.5, 0.5, 50, 300, 50);
[Mx2, My2] = gradient(FOOD);

figure
subplot(121); quiver(Mx1, My1);axis ij
subplot(122); quiver(Mx2, My2);axis ij

%% 2c
% Your answer here:
% The commands for visulization are shown below. For publish purposes and
% to save time, I comment them out.
% Both low: the coli tends to run for a long time before tumbling and
% tumble for a long time before running. 
% Both high: the coli tends to run for a short time before tumbling and
% tumble for a short time before running. 
% One high one low: (1)the coli tends to run for a short time before tumbling and
% tumble for a long time before running, or (2) the coli tends to run for a 
% long time before tumbling and tumble for a short time before running. 
% One zero: the coli either runs or tumbles forever.


%Your code here
%Both low
%[HISTM, FOOD] =    acoli_hist([0,0], 1, 0.1, 0.1, 50, 300, 50,1);
%Both high
%[HISTM, FOOD] =    acoli_hist([0,0], 1, 0.9, 0.9, 50, 300, 50,1);
%rttp high, ttrp low 
%[HISTM, FOOD] =    acoli_hist([0,0], 1, 0.9, 0.1, 50, 300, 50,1);
%rttp low, ttrp high  
%[HISTM, FOOD] =    acoli_hist([0,0], 1, 0.1, 0.9, 50, 300, 50,1);
%rttp 0
%[HISTM, FOOD] =    acoli_hist([0,0], 1, 0, 0.5, 50, 300, 50,1);
%ttrp 0 
%[HISTM, FOOD] =    acoli_hist([0,0], 1, 0.5, 0, 50, 300, 50,1);





%% 2d
% Your answer here:
% By observing the visualizations above, in these following cases the coli 
% tend to move around in the high concentrated areas. These cases include:
% both high, rttp high ttrp low, and when ttrp 0. 
% The reason for this is that the coli needs to tumble to find the good
% concentration of food, and running too much does not quite help.


%Your code here
%Both high
[HISTM, FOOD] =    acoli_hist([0,0], 1, 0.9, 0.9, 50, 300, 50);
figure
subplot(131); imagesc(HISTM);
subplot(132); contour(HISTM);axis ij;
subplot(133); surf(HISTM);axis ij;

%rttp high, ttrp low 
[HISTM, FOOD] =    acoli_hist([0,0], 1, 0.9, 0.1, 50, 300, 50);
figure
subplot(131); imagesc(HISTM);
subplot(132); contour(HISTM);axis ij;
subplot(133); surf(HISTM);axis ij;

%ttrp 0 
[HISTM, FOOD] =    acoli_hist([0,0], 1, 0.5, 0, 50, 300, 50);
figure
subplot(131); imagesc(HISTM);
subplot(132); contour(HISTM);axis ij;
subplot(133); surf(HISTM);axis ij;

%% 2e
% Your answer here: 
% The three cases are discussed below. 
% I simply ignored any steps that occur outside this area for the average.
% The average value of the food function seen by A.coli for each of the
% cases are printed out, and it turned out that the cases with low or zero
% ttrp and high rttp performed the best. 
% The reason for this is that with low ttrp, the coli tends to stick around
% the area of greater function value once the area is reached. 

%Your code here

%Both high
[HISTM, FOOD] =    acoli_hist([0,0], 2, 0.9, 0.9, 50, 300, 50);
figure
subplot(131); imagesc(HISTM);
subplot(132); contour(HISTM);axis ij;
subplot(133); surf(HISTM);axis ij;
average1 = sum(sum(HISTM .* FOOD)); 
disp(average1);

%rttp high ttrp low
[HISTM, FOOD] =    acoli_hist([0,0], 2, 0.9, 0.1, 50, 300, 50);
figure
subplot(131); imagesc(HISTM);
subplot(132); contour(HISTM);axis ij;
subplot(133); surf(HISTM);axis ij;
average2 = sum(sum(HISTM .* FOOD)); 
disp(average2);


%ttrp 0
[HISTM, FOOD] =    acoli_hist([0,0], 2, 0.5, 0, 50, 300, 50);
figure
subplot(131); imagesc(HISTM);
subplot(132); contour(HISTM);axis ij;
subplot(133); surf(HISTM);axis ij;
average3 = sum(sum(HISTM .* FOOD)); 
disp(average3);

%Old value
[HISTM, FOOD] =    acoli_hist([0,0], 2, 0.5, 0.5, 50, 300, 50);
figure
subplot(131); imagesc(HISTM);
subplot(132); contour(HISTM);axis ij;
subplot(133); surf(HISTM);axis ij;
average4 = sum(sum(HISTM .* FOOD)); 
disp(average1);




%% 2f
% Your answer here
% Different cases are run below and plots are generated. More bins is not
% necessarily better because if bin numbers are too large, the areas
% covered by a square element will be small and it will be harder to spot a
% central tendency. 

%Your code here
[HISTM1, FOOD] =    acoli_hist([0,0], 1, 0.5, 0.5, 50, 300, 10);

[HISTM2, FOOD] =    acoli_hist([0,0], 1, 0.5, 0.5, 50, 300, 50);

[HISTM3, FOOD] =    acoli_hist([0,0], 1, 0.5, 0.5, 50, 300, 100);

[HISTM4, FOOD] =    acoli_hist([0,0], 1, 0.5, 0.5, 50, 300, 500);
figure
subplot(141); imagesc(HISTM1);axis xy;
subplot(142); imagesc(HISTM2);axis xy;
subplot(143); imagesc(HISTM3);axis xy;
subplot(144); imagesc(HISTM4);axis xy;

%% 2g
% Your answer here
% If the Paolina image is selected, the food will be represented by the
% numbers in the matrix im, where higher numbers get darker when
% visualized. Say if we still choose the initial position to be (0,0), and
% the parameters to be very high rttp and very low ttrp, then I expect the
% coli to move forward in the image, make a few tumbles, and then always
% tumbles around the region once the coli reaches the high concentration
% parts, e.g. the border of the face. 

%%
close all