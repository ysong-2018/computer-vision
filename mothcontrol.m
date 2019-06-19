function control_function = mothcontrol

  % Any setup code goes here. In particular, this is where you
  % will create your weight matrix.

  weight = ones(256);
  %disp(weight);
  
  

  
  % FIRST PART OF YOUR CODE HERE


  % The next line outputs a handle to firing_function below.
  control_function = @firing_function;
  
  % (If you don't know what function handles are, they are basically a way of
  % storing functions in variables, so that a function can be called without 
  % needing to know its name -- similar to function pointers if you know C.)
  

  function firing_pattern = firing_function(image)

    % Given the image (a 2-D matrix), evaluate the network response,
    % i.e., apply the matrix you created above to the input image. 
    % Your output is given by the "firing_pattern" output variable, 
    % which should be a 4x1 boolean vector (see the help for mothsim
    % for information about which entries correspond to which directions).
    %disp(image);

    image = image(:); %vectorize it first
    %disp(image);


    % SECOND PART OF YOUR CODE - it should consist of two parts 
    % (possibly only two lines of code):
    % 1) multiply your matrix from above by the input image vector
    % 2) apply the threshold, storing the result into the firing_pattern vector
    
    weight_up = ones(256);
    weight_down = ones(256);
    weight_right = ones(256);
    weight_left = ones(256);
  
    weight_up(129:256, 1:256 ) = -1;
    weight_down(1:128, 1:256 ) = -1;
    weight_right(1:256, 1:128 ) = -1;
    weight_left(1:256, 129:256 ) = -1;

  
  %disp(weight_up);
  
  
    weight_up = weight_up(:);
    weight_down = weight_down(:);
    weight_right = weight_right(:);
    weight_left = weight_left(:);

    up = (sum(weight_up.*image));
    down = (sum(weight_down.*image));
    right = (sum(weight_right.*image));
    left = (sum(weight_left.*image));
    
    %disp(up,down,right,left);
    
    up = up>0;
    down = down>0;
    right = right>0;
    left = left>0;
    
    disp(up);
    disp(down);
    disp(right);
    disp(left);
    disp("--------");
    

    firing_pattern = [up;down;right;left]; %(replace this line with your own code)

    %firing_pattern = [false;false;false;false]; %(replace this line with your own code)


    
  end
end
