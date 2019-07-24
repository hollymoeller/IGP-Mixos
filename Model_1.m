function [tset,bio] = Model_1(t_end,ic,pmax1,pmax2,k1,k2,H1,H2,I_in,L1,L2,A,B)
%This is the model of HV Moeller and MG Neubert
%Based on Huisman and Weissman (1994) Ecology
%Model written 3 March 2015
%Initially coded 3 March 2015

%This model describes the competitive interactions between two species
%growing in a water column of depth z.
%These species compete for light.

%We have added intraguild predation:
%Species 1 can be consumed by Species 2, with an attack rate A, and
%conversion efficiency B/A.

tset = linspace(0,t_end,1000);

%[~,bio] = ode45(@bioldynam,tset,ic);

options = odeset('RelTol',1E-9,'AbsTol',1E-9);
[~,bio] = ode23s(@bioldynam,tset,ic,options);

    function [model] = bioldynam(t,y)
       model = zeros(2,1);
       kappa = k1*y(1) + k2*y(2);
       %First equation: Species 1 (prey)'s dynamics
       model(1) = pmax1*y(1)/kappa*log((H1+I_in)/(H1+I_in*exp(-kappa)))...
           - L1*y(1) - A*y(1)*y(2);
       %Second equation: Species 2 (predator)'s dynamics
       model(2) = pmax2*y(2)/kappa*log((H2+I_in)/(H2+I_in*exp(-kappa)))...
           - L2*y(2) + B*A*y(1)*y(2);
    end


end