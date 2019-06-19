function dispmatrix = disparity(leftimg, rightimg, halfpatchsize, maxdisp)
%DISPARITY Calculates the disparity for a pair of rectified stereo images.
%   DISPMATRIX = DISPARITY(LEFTIMG, RIGHTIMG, HALFPATCHSIZE, MAXDISP) calculates
%   the pointwise disparity between a pair of rectified stereo images.
%
%   LEFTIMG is the left image of the stereo pair, while RIGHTIMG is the
%   right image.
%
%   HALFPATCHSIZE specifies half of the patch size (M in the problem
%   statement in the homework). The patch size actually used ends up being
%   2*HALFPATCHSIZE + 1.
%
%   MAXDISP is an optional parameter specifying the maximum absolute
%   disparity to be tested for. If omitted, it defaults to HALFPATCHSIZE.
%
%   The output argument DISP should contain the disparity between the two
%   images at every point where a valid patch comparison can be made. Due
%   to boundary effects, the size of DISPMATRIX will be
%
%       size(leftimg) - 2*halfpatchsize - [0, 2*maxdisp]

if ~exist('maxdisp', 'var')
  maxdisp = halfpatchsize;
end

%get image dimensions
nrows = size(leftimg, 1);
ncols = size(leftimg, 2);

if any([nrows, ncols] ~= size(rightimg))
  error('Left and right images aren''t of the same size');
end

%make sure you understand why the disparity matrix will be smaller than the input images!
dispmatrix = zeros(size(leftimg) - 2*halfpatchsize - [0, 2*maxdisp]);

%range of values for the central patch positions
y0s = (halfpatchsize + 1):(nrows - halfpatchsize);
x0s = (halfpatchsize + maxdisp + 1):(ncols - halfpatchsize - maxdisp);




% this loops through all the possibilities of x0s and y0s
for j=1:length(x0s)
    for i=1:length(y0s)

        %get actual (x0,y0) position in the image
        x0 = x0s(j);
        y0 = y0s(i);

        maxdx = -maxdisp;
        maxsim = -1;
        for dx = -maxdisp:maxdisp
            leftslice = leftimg(-halfpatchsize + y0 : halfpatchsize+y0,  - halfpatchsize + x0: halfpatchsize + x0);
            rightslice = rightimg(  -halfpatchsize + y0:halfpatchsize + y0,  dx - halfpatchsize +x0: dx + halfpatchsize +x0);
            dot_prod = sum(sum(leftslice.*rightslice));
            norml = sqrt(sum(sum(leftslice.*leftslice)));
            normr = sqrt(sum(sum(rightslice.*rightslice)));
            similarity = dot_prod / (norml * normr);
            if similarity > maxsim
                maxsim = similarity;
                maxdx = dx;
            end
   
        end
        dispmatrix(i,j) = maxdx; %disparity value found for current position;
    end
end

end