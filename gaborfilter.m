function gb = gaborfilter(theta, wavelength, phase, sigma, aspect, ksize)
%GABORFILTER Create a Gabor filter.
%   GB = GABORFILTER(THETA, WAVELENGTH, PHASE, SIGMA, ASPECT, KSIZE)
%   creates a Gabor filter GB with orientation THETA (in radians),
%   wavelength WAVELENGTH (in pixels), phase offset PHASE (in radians),
%   envelope standard deviation SIGMA, aspect ratio ASPECT, and dimensions
%   KSIZE x KSIZE. KSIZE is an optional parameter, and if omitted default
%   dimensions are selected.

if ~exist('ksize', 'var')
  ksize = 8*sigma*aspect;
end

if numel(ksize) == 1
  ksize = [ksize ksize];
end

xmax = floor(ksize(2)/2);
xmin = -xmax;
ymax = floor(ksize(1)/2);
ymin = -ymax;

[xs, ys] = meshgrid(xmin:xmax, ymax:-1:ymin);

% Your code here

%ˆx = xcosθ + ysinθ and ˆy = −xsinθ + ycosθ

x = xs*cos(theta) + ys*sin(theta);
y = -xs*sin(theta) + ys*cos(theta);

G = sin((2*pi/wavelength)*y + phase).*exp(-(x.^2/aspect^2 + y.^2)/(2*sigma));
G1 = G - mean(mean(G));
%disp(G1);
G2 = G1 / sqrt(sum(sum(G1.^2)));
%disp(G2);
disp(sum(sum(G2.^2)));
disp(mean(mean(G2)));
gb = G2;





