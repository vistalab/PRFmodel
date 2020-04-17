% NoiseCircleTest (Simulation:  Take a circle and make two noisy estimates of
% its diameter.  (D1 + noise)/(D2 + noise).  This ratio will always be centered
% on 1.  But now, use the same data to estimate
% (max(D1+noise,D2+noise))/min((D1+noise,D2+noise).  This will always be > 1.
% How much greater?  If your estimate is < (1.2?  1.3?) then the data are
% consistent with a circle.



%% Noisy Circle Tests.m
% We have some results in mrVista with simulated data that are not perfect but
% that can be explained by the very nature of the calculation. 




    radius      = 2;  % degrees
    noiseFactor = 0.7;
    % rng(44444,'twister')
    n1          = noiseFactor * randn(1000,1);
    % rng(44322,'twister')
    n2          = noiseFactor * randn(1000,1);

    R1 = radius + n1;
    R2 = radius + n2;

    freeDiff  = R1 - R2;
    freeRatio = R1 ./ R2;
    % mean(freeRatio)
    % std(freeRatio)

limRatio = max([R1,R2]') ./ min([R1,R2]');
mean(limRatio)
% std(limRatio)
median(limRatio)




% Why is the lower part like this? When the noise is similar to the radio
radius       = 0.5;  % degrees
noiseFactor = 0.7;
% rng(44444,'twister')
n1          = noiseFactor * randn(1000,1);
% rng(44322,'twister')
n2          = noiseFactor * randn(1000,1);
R1 = max(0.1,radius + n1);
    R2 = max(0.1,radius + n2);

    freeDiff  = R1 - R2;
    freeRatio = R1 ./ R2;
    % mean(freeRatio)
    % std(freeRatio)

limRatio = max([R1,R2]) ./ min([R1,R2]);
mean(limRatio)
% std(limRatio)
median(limRatio)

