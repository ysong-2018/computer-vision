%% Homework 3
% * Name: Audrey Song
% * NetId: ys585

%IMPORTANT: Please do not modify any of the question headings!
%Also, please do not modify any of the file names.

%% 1a (optional)
% Your answer here
% This is a linear combination because it is is an expression constructed 
% from a set of terms by multiplying each term by a constant and adding 
% the results. 

%Your code here
u= [1 1];
v = [3 -1];
w = 2*u+1.5*v;

%% 1b (optional)
% Your answer here
% For this one, calculate by hand, alpha'=7, beta'=-1

%Your code here

%% 1c (optional)
% Your answer here
% As shown in code

%Your code here
m = transpose([u;v]);
result = inv(m)*transpose([4 8]);

%% 1d (optional)
% Your answer here
% They are independent because one is not scalar multiple of the other

%Your code here

%% 1e (optional)
% Your answer here
% When calculated by hand, there is no valid solution

%Your code here


%% 1f (optional)
% Your answer here
% They are not orthogonal because not all of the pairs have dot product of
% zero
% They form a basis because the vectors are mutually independent and thus
% form a spanning set in R3

%Your code here
x =[1 0 0];
y = [1 -1 -1];
z = [0 3 -1];

%% 1g (optional)
% Your answer here
% They form a basis for R3 and it is an orthogonal basis. Let the column 
% vectors be i1, i2, i3, then the linear combination is 2*i1+i2+7*i3. The
% scalars we multiply are just the entries of the vector

%Your code here

%% 1h (optional)
% Your answer here
% They do form a basis for R2. 
% An orthonormal basis can be acheived by doing matrix reductions and
% giving vectors [1 0] and [0 1]. 

%Your code here

%% 1i (optional)
% Your answer here
% alpha=-1, beta=3

%Your code here

%%
close all

%% 2a (optional)
% Your answer here

%Your code here
t = linspace(0,2*pi,100);
y1 = sin(t);
plot(t,y1)

%% 2b (optional)
% Your answer here

%Your code here
y2 = sin(2*t);
plot(t,y2);hold on;



%% 2c (optional)
% Your answer here

%Your code here

%% 2d (optional)
% Your answer here

%Your code here

%% 2e (optional)
% Your answer here

%Your code here

%% 2f (optional)
% Your answer here

%Your code here

%% 2g (optional)
% Your answer here

%Your code here

%% 2h (optional)
% Your answer here

%Your code here

%% 2i (optional)
% Your answer here

%Your code here

%%
close all

%% 3a
% Your answer here
% Each manually calculated convolutions are put into arrays con1, conv2,
% and conv3 respectively. The stem plots are shown below.

%Your code here
conv1 = [16    32    48    64    80    64    48    32    16];
conv2 = [1     3     5     8    12    14    11     5     1];
conv3 = [12    20    12     4     0    12    20    12     4];

figure
subplot(131);stem(conv1);axis equal; grid on;
subplot(132);stem(conv2);axis equal; grid on;
subplot(133);stem(conv3);axis equal; grid on;

%%
close all

%% 4a
% Your answer here
% The first part of the below demonstrates the Fast Fourier Transform (FFT)
% of the convolution (conv) of f and g. The second part demonstrates the 
% elementwise product between the FFT of f and the FFT of g. 
% Comparison: The real parts and the imaginary parts of both cases resemble
% so much that difference cannot be spotted by eye. Therefore, the
% assumption that the Fourier transform of the convolution of two
% functions is the product of the Fourier transforms of the functions holds
% under Matlab proof.

%Your code here
x = 0:.5:10;
f = x.*(10-x);
g = x.*(10-x).*(5-x);

h = conv(f,g);
y=fft(h,50);
R = real(y);
I = imag(y);
figure
subplot(121); plot(R);
subplot(122); plot(I);

%%
fft_f = fft(f,50);
fft_g = fft(g,50);
prod_fg = fft_f.* fft_g;
prod_fg_R = real(prod_fg);
prod_fg_I = imag(prod_fg);

figure
subplot(121); plot(prod_fg_R);
subplot(122); plot(prod_fg_I);


%%
close all

%% 5a
% Your answer here
% The filter is created as in blurfilter.m. 
% Why is the normaliztion step desired: when convolution is applied to each
% single subgrid of the matrix, the products of the filter and the subgrids
% will be calculated. If the magnitudes of elements in the filter are too
% large or too small, the product matrix elements will have very large or
% very small magnitutes, making it bad for visualization as gives very
% black or white results. 

%Your code in the blurfilter.m file. Please do NOT rename that file.

%% 5b
% Your answer here
% Filters are created and fouriers transforms are applied. The fourier
% transform of the filter looks like a gaussian function. 

%Your code here
filt = blurfilter(1);
f1 = fft2(filt);
f2 = fft2(filt,20,20);

figure
subplot(131); surf(abs(f1), 'EdgeColor', 'none');
subplot(132); surf(abs(f2), 'EdgeColor', 'none');
subplot(133); surf(fftshift(abs(f2)), 'EdgeColor', 'none');




%% 5c
% Your answer here
% The first part gives a filtered Paolina with a blur filter. The second
% part shows the inverse Fourier transform (ifft2) of the pointwise product
% between the Fourier transforms of Paolina and my filter. The results
% shown gives that these two methods give the same blurring effect. In
% fact, what is done here is exactly same as what is demonstrated in
% problem 4, except that this one takes a form of the equation in 4 such
% that both sides are taken an inverse fourier transform. 

%Your code here
im = imread('Paolina.tiff');
im = im2double(im);

filt = blurfilter(3);
imfiltered = filter2(filt,im);
%surf(imfiltered);

f3 = fft2(filt,520,520);
f4 = fft2(im,520,520);
inversed = ifft2(f3.*f4);

figure
subplot(131); imshow(im);
subplot(132); imshow(imfiltered);
subplot(133); imshow(inversed);



%% 5d
% Your answer here
% As shown from the images, this filter returns mostly the contour lines of
% images and ignores the very sharp or very insignificant changes. 
% As we can see from the 2-dimensional Fourier transform of this filter, we
% can see that frequencies that are too low or too high get rejected.
% Therefore, it is a band pass filter as from the definition.

%Your code here
filt = blurfilter(1,13);
filt1 = filt;
filt = blurfilter(2,13);
filt2 = filt;
filt_diff = filt1 - filt2;
disp(filt_diff);

A(1:13,1:13) = mean(mean(filt_diff));
filt_diff = filt_diff - A;
%disp(filt_diff);

im1 = imread('Loki.tiff');
im1 = im2double(im1);

im2 = imread('Steve.tiff');
im2 = im2double(im2);

filtered1 = filter2(filt_diff,im1);
filtered2 = filter2(filt_diff,im2);

figure
subplot(121); imshow(filtered1,[]);
subplot(122); imshow(filtered2,[]);

%%
ft = fft2(filt_diff,15,15);
figure
subplot(121); surf(abs(ft), 'EdgeColor', 'none');
subplot(122); surf(fftshift(abs(f2)), 'EdgeColor', 'none');

%% 5e
% Your answer here
% The unsharpmask.m is implemented. 
% Different cases of pictures are created below and visulized. The
% unsharpmask filter is sharpening the edges of the image. This is because
% the mask identifies the sharpening corners and adding a scaled mask back
% to the original image intensifies the edges. 

%Write your function using the unsharpmask.m file. Please do NOT rename that file.

%Your code here
im = imread('Paolina.tiff');
im = im2double(im);
im1 = imread('Loki.tiff');
im1 = im2double(im1);
im2 = imread('Steve.tiff');
im2 = im2double(im2);

sigma = 1;
amount = 1;

%%
filtered = unsharpmask(im, sigma, amount);
subplot(121); imshow(im);
subplot(122); imshow(filtered);
%%
filtered = unsharpmask(im1, sigma, amount);
subplot(121); imshow(im1);
subplot(122); imshow(filtered);
%%
filtered = unsharpmask(im2, sigma, amount);
subplot(121); imshow(im2);
subplot(122); imshow(filtered);

%% 5f
% Your answer here
% It is possible to design an impulse response you could convolve against 
% the image to perform unsharp masking in a single step. 
% final_im = im + const(im - im*gauss)
% = im + const im - const im*gauss
% im = im*delta
% im const = im*(const delta)
% Plug back in, we get 
% final_im = im*(delta + const delta - const delta*gauss)
% All we need to do is to implement the part given in the parenthesis

%Your code here (if any)

filt = blurfilter(1,13);
gauss = filt;
delta = zeros(13);
delta(7,7)=1;
amount = 1;
conv_in = conv(delta,gauss);
inside = delta + amount*delta + amount*conv_in;
result = conv(im,inside);



%%
close all


