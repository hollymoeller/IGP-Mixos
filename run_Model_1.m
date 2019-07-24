% 3 March 2015
% Coded by HV Moeller

% This function specifies parameters for and runs "Model_1", a
% Huisman-Weissman style model for two competing species of phytoplankton
% growing in a water column. We (Neubert & Moeller) have introduced 
% intraguild predation, in which species 2 can eat (and grows
% from eating) species 2

% PARAMETERS:
I_in = 50;        %Input light to the system
pmax1 = 1; %maximum photosynthesis rate of species 1
pmax2 = .9; %maximum photosynthesis rate of species 2
k1 = .01; %per-cell absorptivity of species 1
k2 = .015; %per-cell absorptivity of species 2
H1 = 1; %half-saturation irradiance of species 1
H2 = 5; %half-saturating irradiance of species 2
L1 = .1; %carbon loss rate of species 1
L2 = .5; %carbon loss rate of species 2
A = 0.003; %attack rate of species 1 on species 2
B = 0.0005; %rescaling, to convert species 1 biomass to species 2 biomass
t_end = 300; %simulation time
ic = [100,300]; %initial conditions, [amount of prey, amount of predator]

% PI CURVES
PI1 = @(I) pmax1 .* (I./(H1+I)) - L1;
PI2 = @(I) pmax2 .* (I./(H2+I)) - L2;

% MODEL CALL:
[tset,bio] = Model_1(t_end,ic,pmax1,pmax2,k1,k2,H1,H2,I_in,L1,L2,A,B); 
I_out = I_in*exp(-(k1*bio(end,1)+k2*bio(end,2)));

% PLOTTING RESULTS
% PLOTTING PARAMETERS
w1col = [0,0,.8];
w2col = [.8,0,0];
coexcol = [.5,0,.5];

% MAKE THE FIGURE 
figure(1)
clf(1)
figure(1)

subplot(2,1,1) %Species Dynamics from this model run
plot(tset,bio(:,1),'Color',w1col)
ylabel('Species 1')
hold on
plot(tset,bio(:,2),'Color',w2col)
ylabel('Species 2')
ylim([0,max(max(bio))]*1.1)
title('Population Dynamics')
xlabel('Time')
ylabel('Population Size')

subplot(2,1,2) %PI curves
plot(linspace(0,I_in*1.5),PI1(linspace(0,I_in*1.5)),'color',w1col)
hold on
plot(linspace(0,I_in*1.5),PI2(linspace(0,I_in*1.5)),'Color',w2col)
yL = get(gca,'YLim');
line([I_out,I_out],yL,'Color','k')
line([I_in,I_in],yL,'Color','k','LineStyle','--')
ylabel('Photosynthetic Rate')
xlabel('Irradiance')
title('P vs. I curves')



