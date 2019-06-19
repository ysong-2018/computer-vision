%% Homework 2
% * Name: Audrey Song
% * NetId: ys585

%IMPORTANT: Please do not modify any of the question headings!

%% 1a
% Your answer here
% For 0nA, the membrane potential gradually decreases with no spike
% For 1nA, the membrane potential remains the same
% For 10nA, the membrane potential has periodic spikes and decreases in
% each period

%Your code here

% 0nA

nsim = nsimulator(1);
ie = 0; % nA

% Set initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, -60);
% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);


% Animate the behavior (display membrane potential vs time)
nsim.animate_integration(nsim, 1, 0.2, 200, [1]);

%%
% 1nA

nsim = nsimulator(1);
ie = 1; % nA

% Set initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, -60);
% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);


% Animate the behavior (display membrane potential vs time)
nsim.animate_integration(nsim, 1, 0.2, 200, [1]);

%%
% 10nA

nsim = nsimulator(1);
ie = 10; % nA

% Set initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, -60);
% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);


% Animate the behavior (display membrane potential vs time)
nsim.animate_integration(nsim, 1, 0.2, 200, [1]);

%% 1b
% Your answer here
% When the initial potential is -80mV: 
% For 0nA, the membrane potential gradually increases to -70 with no spike
% For 1nA, the membrane potential gradually increases to -60 with no spike
% For 10nA, the membrane potential has periodic spikes and increases in
% each period
% When the initial potential is 0mV: 
% For 0nA, the membrane potential gradually increases to -70 with an
% initial spike
% For 1nA, the membrane potential gradually increases to -60 with an
% initial spike
% For 10nA, the membrane potential has periodic spikes and increases in
% each period


%Your code here
% -80mV
% 0nA

nsim = nsimulator(1);
ie = 0; % nA

% Set initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, -80);
% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);


% Animate the behavior (display membrane potential vs time)
nsim.animate_integration(nsim, 1, 0.2, 200, [1]);

%%
% 1nA

nsim = nsimulator(1);
ie = 1; % nA

% Set initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, -80);
% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);


% Animate the behavior (display membrane potential vs time)
nsim.animate_integration(nsim, 1, 0.2, 200, [1]);

%%
% 10nA

nsim = nsimulator(1);
ie = 10; % nA

% Set initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, -80);
% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);


% Animate the behavior (display membrane potential vs time)
nsim.animate_integration(nsim, 1, 0.2, 200, [1]);


%%
% 0mV
% 0nA

nsim = nsimulator(1);
ie = 0; % nA

% Set initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, 0);
% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);


% Animate the behavior (display membrane potential vs time)
nsim.animate_integration(nsim, 1, 0.2, 200, [1]);

%%
% 1nA

nsim = nsimulator(1);
ie = 1; % nA

% Set initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, 0);
% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);


% Animate the behavior (display membrane potential vs time)
nsim.animate_integration(nsim, 1, 0.2, 200, [1]);

%%
% 10nA

nsim = nsimulator(1);
ie = 10; % nA

% Set initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, 0);
% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);


% Animate the behavior (display membrane potential vs time)
nsim.animate_integration(nsim, 1, 0.2, 200, [1]);


%% 1c
% Your answer here
% The different cases are run below: 
% From the outputs, we can see that the membrane potentials of neuron 1 has
% a negative impact on that of neuron 2



%Your code here
% 0nA
nsim = nsimulator(2);
ie = 2.5; % nA

% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);
nsim = nsim.set_Ie(nsim, 2, 0);

% Set the initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, -60);
nsim = nsim.set_mpotential(nsim, 2, -70);

% Specify connections
% nsim = nsim.make_excitatory_connection(nsim, 2, 1);
nsim = nsim.make_inhibitory_connection(nsim, 1, 2);

% Simulate and plot
norecordtime = 0.0; % Integration time (seconds) when not recording 
                    % membrane potentials
recordtime = 0.1; % Integration time for recorded values
nsim = nsim.simulate_timeinterval(nsim, norecordtime, recordtime);

% Display the recorded values
%
% Notice that no current is injected into neuron 3 and that the neuron is
% completely isolated from the other two (no connections).

figid = 2;
nsim.plot_membranepotentials(nsim, figid, [1 2])

%%
% 1nA
nsim = nsimulator(2);
ie = 2.5; % nA

% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);
nsim = nsim.set_Ie(nsim, 2, 1);

% Set the initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, -60);
nsim = nsim.set_mpotential(nsim, 2, -70);

% Specify connections
% nsim = nsim.make_excitatory_connection(nsim, 2, 1);
nsim = nsim.make_inhibitory_connection(nsim, 1, 2);

% Simulate and plot
norecordtime = 0.0; % Integration time (seconds) when not recording 
                    % membrane potentials
recordtime = 0.1; % Integration time for recorded values
nsim = nsim.simulate_timeinterval(nsim, norecordtime, recordtime);

% Display the recorded values
%
% Notice that no current is injected into neuron 3 and that the neuron is
% completely isolated from the other two (no connections).

figid = 2;
nsim.plot_membranepotentials(nsim, figid, [1 2])

%%
% 10nA
nsim = nsimulator(2);
ie = 2.5; % nA

% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);
nsim = nsim.set_Ie(nsim, 2, 10);

% Set the initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, -60);
nsim = nsim.set_mpotential(nsim, 2, -70);

% Specify connections
% nsim = nsim.make_excitatory_connection(nsim, 2, 1);
nsim = nsim.make_inhibitory_connection(nsim, 1, 2);

% Simulate and plot
norecordtime = 0.0; % Integration time (seconds) when not recording 
                    % membrane potentials
recordtime = 0.1; % Integration time for recorded values
nsim = nsim.simulate_timeinterval(nsim, norecordtime, recordtime);

% Display the recorded values
%
% Notice that no current is injected into neuron 3 and that the neuron is
% completely isolated from the other two (no connections).

figid = 2;
nsim.plot_membranepotentials(nsim, figid, [1 2])


%% 1d
% Your answer here
% The three cases are run below. For an excitatory connection, there is not
% much effect of neuron 1 on neuron 2. 


%Your code here
% 0nA
nsim = nsimulator(2);
ie = 2.5; % nA

% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);
nsim = nsim.set_Ie(nsim, 2, 0);

% Set the initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, -60);
nsim = nsim.set_mpotential(nsim, 2, -70);

% Specify connections
nsim = nsim.make_excitatory_connection(nsim, 2, 1);

% Simulate and plot
norecordtime = 0.0; % Integration time (seconds) when not recording 
                    % membrane potentials
recordtime = 0.1; % Integration time for recorded values
nsim = nsim.simulate_timeinterval(nsim, norecordtime, recordtime);

% Display the recorded values
%
% Notice that no current is injected into neuron 3 and that the neuron is
% completely isolated from the other two (no connections).

figid = 2;
nsim.plot_membranepotentials(nsim, figid, [1 2])

%%
% 1nA
nsim = nsimulator(2);
ie = 2.5; % nA

% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);
nsim = nsim.set_Ie(nsim, 2, 1);

% Set the initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, -60);
nsim = nsim.set_mpotential(nsim, 2, -70);

% Specify connections
nsim = nsim.make_excitatory_connection(nsim, 2, 1);

% Simulate and plot
norecordtime = 0.0; % Integration time (seconds) when not recording 
                    % membrane potentials
recordtime = 0.1; % Integration time for recorded values
nsim = nsim.simulate_timeinterval(nsim, norecordtime, recordtime);

% Display the recorded values
%
% Notice that no current is injected into neuron 3 and that the neuron is
% completely isolated from the other two (no connections).

figid = 2;
nsim.plot_membranepotentials(nsim, figid, [1 2])

%%
% 10nA
nsim = nsimulator(2);
ie = 2.5; % nA

% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);
nsim = nsim.set_Ie(nsim, 2, 10);

% Set the initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, -60);
nsim = nsim.set_mpotential(nsim, 2, -70);

% Specify connections
nsim = nsim.make_excitatory_connection(nsim, 2, 1);

% Simulate and plot
norecordtime = 0.0; % Integration time (seconds) when not recording 
                    % membrane potentials
recordtime = 0.1; % Integration time for recorded values
nsim = nsim.simulate_timeinterval(nsim, norecordtime, recordtime);

% Display the recorded values
%
% Notice that no current is injected into neuron 3 and that the neuron is
% completely isolated from the other two (no connections).

figid = 2;
nsim.plot_membranepotentials(nsim, figid, [1 2])

%% 1e
% Your answer here
% Initially, I do not see an offset. After setting the membrane potentials
% randomly, the offset appears. Stimulating several times gives different
% offsets since the initial potentials are all set randomly. 
% For inhibitory connections, neurons tend to synchronize their spikes

%Your code here
nsim = nsimulator(2);
ie = 2.5; % nA

% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);
nsim = nsim.set_Ie(nsim, 2, 2.5);

% Get random potentials between -70 ~ -50
a = -70;
b = -50;
r1 = (b-a).*rand(1,1) + a;
r2 = (b-a).*rand(1,1) + a;

% Set the initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, r1);
nsim = nsim.set_mpotential(nsim, 2, r2);

% Specify connections
nsim = nsim.make_excitatory_connection(nsim, 2, 1);
nsim = nsim.make_excitatory_connection(nsim, 1, 2);


% Simulate and plot
norecordtime = 1.0; % Integration time (seconds) when not recording 
                    % membrane potentials
recordtime = 0.1; % Integration time for recorded values
nsim = nsim.simulate_timeinterval(nsim, norecordtime, recordtime);

% Display the recorded values
%
% Notice that no current is injected into neuron 3 and that the neuron is
% completely isolated from the other two (no connections).

figid = 2;
nsim.plot_membranepotentials(nsim, figid, [1 2])

%%
% 
nsim = nsimulator(2);
ie = 2.5; % nA

% Define injected current
nsim = nsim.set_Ie(nsim, 1, ie);
nsim = nsim.set_Ie(nsim, 2, 2.5);

% Get random potentials between -70 ~ -50
a = -70;
b = -50;
r1 = (b-a).*rand(1,1) + a;
r2 = (b-a).*rand(1,1) + a;

% Set the initial membrane potential
nsim = nsim.set_mpotential(nsim, 1, r1);
nsim = nsim.set_mpotential(nsim, 2, r2);

% Specify connections
nsim = nsim.make_inhibitory_connection(nsim, 2, 1);
nsim = nsim.make_inhibitory_connection(nsim, 1, 2);


% Simulate and plot
norecordtime = 1.0; % Integration time (seconds) when not recording 
                    % membrane potentials
recordtime = 0.1; % Integration time for recorded values
nsim = nsim.simulate_timeinterval(nsim, norecordtime, recordtime);

% Display the recorded values
%
% Notice that no current is injected into neuron 3 and that the neuron is
% completely isolated from the other two (no connections).

figid = 2;
nsim.plot_membranepotentials(nsim, figid, [1 2])


%%
close all

%% 2a
% Your answer here
% Doing the multiplication gives the same vector

%Your code here

vec = [1 2 3 4 5];
A = eye(5);
B = A*transpose(vec);
disp(B);


%% 2b
% Your answer here
% Each element of the vector becomes twice the times the original one


%Your code here

vec = [1 2 3 4 5];
C = 2*eye(5);
D = C*transpose(vec);
disp(D);

%% 2c
% Your answer here
% This time, each element of the vector becomes the sum of the element of
% the input vector plus the neighbors of the input vector

%Your code here
vec = [1 2 3 4 5];

E = eye(5);
for i = 1:4
    E(i,i+1) = 1;
end

for i = 2:5
    E(i,i-1) = 1;
end

F = E*transpose(vec);
disp(F);

%% 2d
% Your answer here
% The resulting matrix is G. This is because the diagonal terms take care
% of the original terms of each vector element and the diagonals above and
% below the main diagonal take care of the neighboring terms

%Your code here

G = eye(10);
for i = 1:9
    G(i,i+1) = -0.2;
end

for i = 2:10
    G(i,i-1) = -0.2;
end

disp(G);

%% 2e
% Your answer here
% Plotting a using plot command gives a function of y=10 where x=1~5 and 
% y=20 where x=6~10
% Plotting the output vector H using plot command gives a zig zag line that
% reflects the impact of the matrix on the original vector
% Visualizing them using imagesc and colormap gray reflects the values of
% each vector on the darkness of the sections in the plot. The results
% match the results of the plot command qualitatively
% The matrix actually sharpens the input vector by enhancing the edges. 
% Consider that when we subtract the values from the left and right, we are 
% in a way exaggerating the differences between adjacent values, which 
% sharpens the image


%Your code here
a= [10 10 10 10 10 20 20 20 20 20];

H = G*transpose(a);


figure
subplot(121); plot(a);
subplot(122); plot(H);

figure
subplot(121); imagesc(a);
subplot(122); imagesc(H(2:9));
colormap gray;



%% 2f
% Your answer here
% I built the matrix J reflecting this feedforward model. 
% This network make each element of the output vector = element of the
% input vector + (1/3)*(sum of neighbors)
% g on here: we are summing the values of adjacent pixels, which in a way 
% acts as an unscaled averaging of the values. This diminishes the apparent
% contrast, thus blurring the image

%Your code here
J = eye(10);
for i = 1:9
    J(i,i+1) = 1/3;
end

for i = 2:10
    J(i,i-1) = 1/3;
end

b = [10 10 10 10 10 20 20 20 20 20];
c = [10 8 10 8 10 8 10 8 10 8];
b2 = J*transpose(b);
c2 = J*transpose(c);


figure
subplot(121); plot(b2);
subplot(122); plot(c2);
colormap gray;




%% 2g
% Your answer here
% To apply the two layers, we apply G first and them apply J to the
% resulting vector from the previous step. If we apply this to vector a, it
% will look like a2 = J*G*transpose(a); 
% The result will give each element a function of summing a certain amount
% of the original vector and also involving two neighboring terms on each
% side. 


%Your code here
d = [10 10 10 10 10 10 10 10 10 10];
K = J*G;
d2 = K*transpose(d);
disp(K);
disp(d2);



%%
close all

%% 3a
% Your answer here
% The input nodes are the pixels of the image from our model moth eye, 
% and whose four output nodes control the activity of the moth's motor 
% neurons so as to make it fly towards light. In my design, I will let the
% moth judge from its eye the image. If the sum of the upper half matrix is
% greater than that of the lower half, then it should move up. If the sum 
% of the upper half matrix is less greater than that of the lower half, 
% then it should move down. If the sum of the right half matrix is
% greater than that of the left half, then it should move right. If the sum
% of the left half matrix is greater than that of the right half, then it 
% should move left. The up, down, left, right boolean values form a vector
% that can be used. 
% Suppose the sum of one half of the matrix is a and the sum of the other
% half of the matrix is b, then for each of the four cases, we apply the
% threshold a>b which is a boolean value. 

%Your code here

%3b - in the mothcontrol.m file.