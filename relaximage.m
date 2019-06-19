function relaxed = relaximage(image, niters)

[x,y] = size(image);

%Preprocessing of image adding a pad to the surroundings
prev = zeros(x+2,y+2);
prev(2:x+1,2:y+1) = image;
prev(1,1) = prev(2,2);
prev(1,y+2) = prev(2,y+1);
prev(x+2,1) = prev(x+1,2);
prev(x+2,y+2) = prev(x+1,y+1);
prev(1,2:y+1) = prev(2,2:y+1);
prev(x+2,2:y+1) = prev(x+1,2:y+1);
prev(2:x+1,1) = prev(2:x+1,2);
prev(2:x+1,y+2) = prev(2:x+1,y+1);

%Let background be 1, object be 0, these are the two lambdas
%We will only calculate the case for lambda=1 since the other will just be
%one minus this

cur = zeros(x+2,y+2);
for n=1:niters
    for a=2:x+1
        for b=2:y+1
            %let point i be our current point(a,b)
            %let the points j be all the eight neighboring points

            %support_l1 = sum_N(r_ij(1,0)*p_j(0) + r_ij(1,1)*p_j(1))
            %let p_j(1) be the pixel value of the image
            %r_ij(1,0)=-1, r_ij(1,1)=1
            %p_j(0) = 1-p_j(1)
            %support_l1 =sum_N(-1*(1-p_j(1)) + 1*p_j(1)) = sum_N(2*p_j(1)-1)

            support = (2*prev(a-1,b-1)-1)+  (2*prev(a-1,b)-1)+  (2*prev(a-1,b+1)-1)+  (2*prev(a,b-1)-1)+  (2*prev(a,b+1)-1)+  (2*prev(a+1,b-1)-1)+  (2*prev(a+1,b)-1)+  (2*prev(a+1,b+1)-1);
            support_scale = support/8;
            trun = prev(a,b)+support_scale;
            if trun<0
                trun=0;
            end
            if trun>1
                trun=1;
            end
            cur(a,b)=trun;
        end
    end

    cur(1,1) = cur(2,2);
    cur(1,y+2) = cur(2,y+1);
    cur(x+2,1) = cur(x+1,2);
    cur(x+2,y+2) = cur(x+1,y+1);
    cur(1,2:y+1) = cur(2,2:y+1);
    cur(x+2,2:y+1) = cur(x+1,2:y+1);
    cur(2:x+1,1) = cur(2:x+1,2);
    cur(2:x+1,y+2) = cur(2:x+1,y+1);

    prev = cur;
end


relaxed = cur;