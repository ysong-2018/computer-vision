%% Homework 5
% * Name: Audrey Song
% * NetId: ys585

%IMPORTANT: Please do not modify any of the question headings!
%Also, please do not modify any of the file names (including this one).

%% 1a
% Your answer here
% The derivative is highest when at x=0, and lowest on both sides. 
% Since the derivative is largest in the middle, where the sigmoid changes
% the fastest, and the derivative on the sides are smallest, where the
% sigmoid changes the slowest: it can be useful when looking at edges.
% Edges are points where image values change the most rapidly, and using
% the sigmoid as well as the derivative, we can easily locate the changing
% point. 

%Your code here
n = 101;
x=linspace(-5,5,n);
sigmoid = sigmf(x,[2 0]);
sig_div = gradient(sigmoid);
plot(x,sigmoid);
hold on;
yyaxis right;
plot(x,sig_div);


%% 1b
% Your answer here
% The instructions are followed and the lines are plotted. 
% It overlaps with the derivative from 1a
% This can be seen as a derivative operation because  [.5 0 -.5] indicates
% the direction of transition, and D is to be applied to the convolution. 

%Your code here
filter = [.5 0 -.5];
sig_grad = conv(sigmoid,filter,'valid');
sig_grad = [0 sig_grad 0];

plot(x,sig_div);
hold on;
yyaxis right;
plot(x,sig_grad);


%% 1c 
% Your answer here
% The output and the original derivative are plotted.

%Your code here
temp = zeros(1,n);

for i=2:n-1 
    %if it is greater than or equal to its left neighbor and greater than
    %its right neighbor, output a 1. else, output a 0.
    if (sig_grad(i) >= sig_grad(i-1)) & (sig_grad(i) > sig_grad(i+1))
        temp(i) = 1;
    else
        temp(i) = 0;
    end 
end

plot(x,sig_grad);
hold on;
yyaxis right;
plot(x,temp);

%% 1d 
% Your answer here
% Follow the given code and plot the answer. Exactly the same as 1c.

%Your code here
xs = 2:n-1; 
xpeaks = sig_grad(1,xs-1) <= sig_grad(1,xs) & sig_grad(1,xs) > sig_grad(1,xs+1);
xpeaks = [0 xpeaks 0];

plot(x,sig_grad);
hold on;
yyaxis right;
plot(x,xpeaks);

%% 1e
% Your answer here
% Different functions are plotted. What is special is that when the second
% derivative crosses zero, the x and y are zero. This is because when x=0,
% the first derivative reaches the maximum, and the derivative of the first
% derivative at this point should be zero. 

%Your code here
sig_div_2 = gradient(gradient(sigmoid));

plot(x,sig_div);
hold on;
yyaxis right;
plot(x,sig_div_2);


%% 1f
% Your answer here
% The similar vector is obtained here. Here we need to use the fact that if
% a function reaches the maximum, the derivative at the maximum point
% equals 0 and derivative around the maximum point decreases. This is how I
% reason my algorithm. 
% However, when trying to use the condition sig_div_2(i)==0, I discovered
% that the gradient function does not give exact zeros, so i changed the
% method to trying to find i where the second derivative around it changes
% the most rapidly.

%Your code here
vec = zeros(1,n);

maxdiff =0;
idx=0;
for i=2:n-1 
    if (sig_div_2(i) <= sig_div_2(i-1)) & (sig_div_2(i) > sig_div_2(i+1)) 
        if  (sig_div_2(i-1) - sig_div_2(i+1)) >maxdiff
            maxdiff = (sig_div_2(i-1) - sig_div_2(i+1));
            idx = i;
        end 
    end 
end

vec(idx) = 1;

plot(x,sig_div_2);
hold on;
yyaxis right;
plot(x,vec);

%%
close all; clear all

%% 2a 
% Your answer here
% The results are shown in the image. 
% From these two sobel operators, we get the two "partial derivatives", by
% seeing in the two separate images the different changes in different
% directions. 

%Your code here
im = imread('Paolina.tiff');
im = im2double(im);

op1 = zeros(3,3);
op2 = zeros(3,3);

op1(1,1) = -1;
op1(2,1) = -2;
op1(3,1) = -1;
op1(1,3) = 1;
op1(2,3) = 2;
op1(3,3) = 1;

op2(1,1) = -1;
op2(1,2) = -2;
op2(1,3) = -1;
op2(3,1) = 1;
op2(3,2) = 2;
op2(3,3) = 1;

par1 = conv2(im,op1,'same');
par2 = conv2(im,op2,'same');

figure
subplot(121);imshow(par1,[]);
subplot(122);imshow(par2,[]);


%% 2b 
% Your answer here
% Code is implemented. By using the magnitude of the gradient, we combine
% the effects of the two partial derivatives. By applying the threshold, we
% only keep high intensity for the edges. 

%Your code here

mag_img = sqrt((par1.^2)+(par2.^2));
thres(1:480,1:512) = 2*mean(mean(mag_img));
clean_mag_img = mag_img > thres ;
imshow(clean_mag_img,[]);

%% 2c 
% Your answer here
% Code is implemented. I can see that the edges get sharpened a lot. the
% purpose of this step is to give a cleaner view of the edge. Compared to
% what is done in 1cd, this problem has the similar idea. The logic is to
% only keep the points where the rate of change acheives maximum, thereby
% sharpening the edges to the most precise. 

%Your code here
[m n] = size(clean_mag_img);
xs = 2:n-1;
ys = 2:m-1;
xpeaks = clean_mag_img(ys,xs-1) <= clean_mag_img(ys,xs) & clean_mag_img(ys,xs) > clean_mag_img(ys,xs+1);
ypeaks = clean_mag_img(ys-1,xs) <= clean_mag_img(ys,xs) & clean_mag_img(ys,xs) > clean_mag_img(ys+1,xs);
clean_mag_img(ys,xs)= clean_mag_img(ys,xs) & ( xpeaks | ypeaks );
imshow(clean_mag_img,[]);

%% 2d
% Your answer here
% The laplacian of the Gaussian can be computed directly and sampled to the
% desired accuracy.
% The results are a little worse than my own sobel detector. The new
% results are a little unclean, having some undesired lines and dots around
% the real edges. 
% The point of applying the gaussian to the laplacian is to give some
% smoothing features. 

%Your code here
edges = edge(im,'log');
imshow(edges,[]);

%%
close all; clear all

%% 3a 
% Your answer here
% Here we can see the canny threshold works the best

%Your code here
im = imread('Circles.tiff');
im = im2double(im);

edge1 = edge(im,'Sobel');
edge2 = edge(im,'Canny');

figure
subplot(121);imshow(edge1,[]);
subplot(122);imshow(edge2,[]);


%% 3b 
% Your answer here
% In this result we can see the approxcanny works the best, giving the
% clearest outer edge, but still does not work ideally to capture the inner
% ones. 

%Your code here
I = imnoise(im, 'gaussian');
%imshow(I);
edge1 = edge(I,'Sobel');
edge2 = edge(I,'Canny');
edge3 = edge(I,'Prewitt');
edge4 = edge(I,'Roberts');
edge5 = edge(I,'log');
edge6 = edge(I,'zerocross');
edge7 = edge(I,'approxcanny');

figure
subplot(171);imshow(edge1,[]);
subplot(172);imshow(edge2,[]);
subplot(173);imshow(edge3,[]);
subplot(174);imshow(edge4,[]);
subplot(175);imshow(edge5,[]);
subplot(176);imshow(edge6,[]);
subplot(177);imshow(edge7,[]);


%% 3c 
% Your answer here
% In this one we tried a few choices. Canny captured all the edges and the
% others captured the outer circle. This is a much better result and for
% sure the blur improves the result. 

%Your code here
sigma = 3;
filtersize = 10*sigma;
blurred = imfilter(I, fspecial('gaussian', filtersize, sigma), 'replicate');
imshow(blurred);

edge1 = edge(blurred,'Sobel');
edge2 = edge(blurred,'Canny');
edge3 = edge(blurred,'Prewitt');
edge4 = edge(blurred,'Roberts');
edge5 = edge(blurred,'log');
edge6 = edge(blurred,'zerocross');
edge7 = edge(blurred,'approxcanny');

figure
subplot(171);imshow(edge1,[]);
subplot(172);imshow(edge2,[]);
subplot(173);imshow(edge3,[]);
subplot(174);imshow(edge4,[]);
subplot(175);imshow(edge5,[]);
subplot(176);imshow(edge6,[]);
subplot(177);imshow(edge7,[]);

%%
close all; clear all

%% 4a
% Your answer here
% The function is implemented and the 16 filters are created as shown
%Please write your function in gaborfilter.m. Please do not rename that file.

%Your code here

% Edge-enhancing filters will have phase offset 0; orientations 0 ~7/4 pi
ef1 = gaborfilter(0, 4, 0, 1, 1.5);
ef2 = gaborfilter(pi/4, 4, 0, 1, 1.5);
ef3 = gaborfilter(pi/2, 4, 0, 1, 1.5);
ef4 = gaborfilter(3*pi/4, 4, 0, 1, 1.5);
ef5 = gaborfilter(pi, 4, 0, 1, 1.5);
ef6 = gaborfilter(5*pi/4, 4, 0, 1, 1.5);
ef7 = gaborfilter(3*pi/2, 4, 0, 1, 1.5);
ef8 = gaborfilter(7*pi/4, 4, 0, 1, 1.5);

%  bright line filters will have phase offset pi/2; orientations 0 ~3/4 pi
blf1 = gaborfilter(0, 4, pi/2, 1, 1.5);
blf2 = gaborfilter(pi/4, 4, pi/2, 1, 1.5);
blf3 = gaborfilter(pi/2, 4, pi/2, 1, 1.5);
blf4 = gaborfilter(3*pi/4, 4, pi/2, 1, 1.5);

%  dark line filters will have phase offset 3*pi/2; orientations 0 ~3/4 pi
dlf1 = gaborfilter(0, 4, 3*pi/2, 1, 1.5);
dlf2 = gaborfilter(pi/4, 4, 3*pi/2, 1, 1.5);
dlf3 = gaborfilter(pi/2, 4, 3*pi/2, 1, 1.5);
dlf4 = gaborfilter(3*pi/4, 4, 3*pi/2, 1, 1.5);
figure
subplot(4,4,1);imagesc(ef1);
subplot(4,4,2);imagesc(ef2);
subplot(4,4,3);imagesc(ef3);
subplot(4,4,4);imagesc(ef4);
subplot(4,4,5);imagesc(ef5);
subplot(4,4,6);imagesc(ef6);
subplot(4,4,7);imagesc(ef7);
subplot(4,4,8);imagesc(ef8);
subplot(4,4,9);imagesc(blf1);
subplot(4,4,10);imagesc(blf2);
subplot(4,4,11);imagesc(blf3);
subplot(4,4,12);imagesc(blf4);
subplot(4,4,13);imagesc(dlf1);
subplot(4,4,14);imagesc(dlf2);
subplot(4,4,15);imagesc(dlf3);
subplot(4,4,16);imagesc(dlf4);






%% 4b
% Your answer here
% Images are shown. Images with the same orientation are located in the
% same columns. Different filters focus on different sides. 

%Your code here

im = imread('Paolina.tiff');
im = im2double(im);

filt1 = conv2(im,ef1);
filt2 = conv2(im,ef2);
filt3 = conv2(im,ef3);
filt4 = conv2(im,ef4);
filt5 = conv2(im,ef5);
filt6 = conv2(im,ef6);
filt7 = conv2(im,ef7);
filt8 = conv2(im,ef8);
filt9 = conv2(im,blf1);
filt10 = conv2(im,blf2);
filt11 = conv2(im,blf3);
filt12 = conv2(im,blf4);
filt13 = conv2(im,dlf1);
filt14 = conv2(im,dlf2);
filt15 = conv2(im,dlf3);
filt16 = conv2(im,dlf4);

figure
subplot(4,4,1);imshow(filt1);
subplot(4,4,2);imshow(filt2);
subplot(4,4,3);imshow(filt3);
subplot(4,4,4);imshow(filt4);
subplot(4,4,5);imshow(filt5);%imshow(ef5,[]);
subplot(4,4,6);imshow(filt6);%imshow(ef6,[]);
subplot(4,4,7);imshow(filt7);%imshow(ef7,[]);
subplot(4,4,8);imshow(filt8);%imshow(ef8,[]);
subplot(4,4,9);imshow(filt9);%imshow(blf1,[]);
subplot(4,4,10);imshow(filt10);%imshow(blf2,[]);
subplot(4,4,11);imshow(filt11);%imshow(blf3,[]);
subplot(4,4,12);imshow(filt12);%imshow(blf4,[]);
subplot(4,4,13);imshow(filt13);%imshow(dlf1,[]);
subplot(4,4,14);imshow(filt14);%imshow(dlf2,[]);
subplot(4,4,15);imshow(filt15);%imshow(dlf3,[]);
subplot(4,4,16);imshow(filt16);%imshow(dlf4,[]);


%% 4c
% Your answer here
% Visual cortex consists of a complex arrangement of cells, with rich 
% interactions within and between layers.
% Gabor filter basically analyzes whether there are any specific frequency 
% content in the image in specific directions in a localized region around 
% the point or region of analysis.
% Therefore, gabors for different orientations can be representative of
% different receptive fields. 
% the computational advantages of the latter representation is that
% different elements of the filters are separated and can be applied in
% more ways. 

%% 4d
% Your answer here
% The frequency domain representation of a Gabor filter is a combination of
% delta functions extending to four directions. Take a random filter, plot
% the fftshift, and get the result. 

%Your code here
surf(fftshift(filt8), 'EdgeColor', 'none');

%%
close all; clear all

%% 5a
% Your answer here
% The disparity matrices give general groups of shades in the images.
% the frontal parallel plane approximation is roughly valid. general trends
% are correctly reflected by the algorithm.
%Please write your function in disparity.m. Please do not rename that file.

%Your code here

left = imread('Scene-Left.tiff');
left = im2double(left);
right = imread('Scene-Right.tiff');
right = im2double(right);
%%
dispmatrix = disparity(left, right, 10);
imagesc(dispmatrix);

% imwrite(mat2gray(dispmatrix),'scene.png');
%%
left = imread('Pentagon-Left.tiff');
left = im2double(left);
right = imread('Pentagon-Right.tiff');
right = im2double(right);
%%
dispmatrix = disparity(left, right, 10);
imagesc(dispmatrix);

% imwrite(mat2gray(dispmatrix),'pentagon.png');

%%
close all; clear all


