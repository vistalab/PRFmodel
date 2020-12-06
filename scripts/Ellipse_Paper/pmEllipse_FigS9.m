function pmEllipse_FigS9
% Figure showing RF spread in a voxel - asking for orientation effect
%
%   And we should turn this into a true function.
%
% Description:
%   Part of the analysis as to whether the expected spread of center
%   positions might produce an oriented pRF.
%
%   TODO:  We need to upgrade to both angular and eccentricity.
%   Specifically, we need to find the magnification factor for the angular
%   dimension and find all the RFs for both eccentricity and angle.
%
%   Also, we should allow that the RF itself have a 2:1 shape.  And we
%   should look in the literature to see whether the support of an oriented
%   RF is circular or elliptical.
%
% See also
%  pmMainEllipseFiguresScript

% x,y single unit RF center positions in degrees
[X,Y] = meshgrid(linspace(-5,5, 200)); 

% RF single unit eccentricity in degrees
E = sqrt(X.^2 + Y.^2); 

% single unit V1 RF radius from macaque (Freeman & Simoncelli, 2011)
R = 0.075*E;

% inverse linear cortical magnification function (deg/mm, Horton and Hoyt, 1991)
invM = @(x) (x+0.75) / 17.3;

% voxel diameter (mm)
d = 2;

% example voxel pRF center (degrees)
x = 2; y = 2; ecc = sqrt(x^2 + y^2);

% inverse cortical magnification at exmaple voxel (deg per mm)
invm = invM(ecc);

% spread of neural RF centers within voxels (deg, voxel radius)
rfSpread = 0.5*d*invm;

% find neurons inside the voxel, ie neurons whose RF centers are within rfSpread of the voxel's center 
inds = find(sqrt((X-x).^2 + (Y-y).^2)  < rfSpread);

% Plot all single unit RFs within voxel
figure(1); clf; 

% do it twice, once zoomed in and once zoomed out
for subs = 1:2
    subplot(1,2,subs); hold on;
    theta = linspace(0,2*pi,100);
    for ii = inds'
        [xi, yi] = pol2cart(theta, R(ii));
        xi = xi + X(ii);
        yi = yi + Y(ii);
        plot(xi,yi, 'k-');
    end

    axis square;
    if subs == 1 % zoom out
        axis(1.5*ecc*[-1 1 -1 1]); 
        ax = gca;
        set(ax, 'FontSize', 12, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')        
    end
end

end

