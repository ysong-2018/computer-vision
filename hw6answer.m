%% Homework 6
% * Name: Audrey Song
% * NetId: ys585

%IMPORTANT: Please do not modify any of the question headings!
%Also, please do not modify any of the file names (including this one).

%% 1
% Your answer here
% Discussion:
% The nodes are each one of the pixels of the images.
% The labels are 0 and 1, representing black and white respectively.
% The initial confidences will be set as the pixel values of the orginal.
% image. (the image is grayscale anyway so we don't need to do anthing on
% it).
% Neighborhood relationships: we choose the surrounding 8 pixels as the set
% of neightborhood nodes.
% Compatibility function: if the labels are the same, then the
% compatibility would be 1. If the labels are different, then the
% compatibility would be -1
% Further detailed discussions are illustrated in the comments of the
% function file. 
% Images:
% Generally, it takes 3-4 iterations for the image to be cleared. However,
% for cactus and parts, where there is a black edge for each one of them,
% it takes quite long for my algorithm to get the black edge converge. This
% is because I added a padding towards the original images for calculations
% purposes, and the padding value takes the closest value from the original
% image.

%Please write your function in relaximage.m. Please do not rename that file.
%Your code here

im1 = imread('cactus.png');
im1 = im2double(im1);

relaxed = relaximage(im1, 105);
imshow(relaxed);

imwrite(mat2gray(relaxed),'cactus105.png');



%%
im2 = imread('parts.png');
im2 = im2double(im2);

relaxed = relaximage(im2, 221);
imshow(relaxed);

imwrite(mat2gray(relaxed),'parts221.png');



%%
im3 = imread('bump.png');
im3 = im2double(im3);

relaxed = relaximage(im3, 10);
imshow(relaxed);

imwrite(mat2gray(relaxed),'bump10.png');


%%
close all; clear all


