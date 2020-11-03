function [area,e] = pmEllipseArea(a,b)
    %This function will calculate the area
    %and the eccentricity of an ellipse if all the
    %conditions for a real ellipse are met.
    %Function usage:
    %pmEllipseArea(a,b)
    %
    %Where a - Semi major axis
    %b - Semi minor axis
    if length(a)==1 && length(b)==1
        if a == 0
            area = pi*a*b;
            e = false;
        elseif a<0 || b<0
            area = false;
            e = false;
        elseif b>a
            area = pi*a*b;
            e = false;
        else
            area = pi*a*b;
            e = sqrt(1-(b^2/a^2));
        end
    end
    
    if length(a) == length(b)
        area = pi * (a .* b);
        e    = sqrt(1-(b .^ 2 ./ a .^ 2));
    else
        error('Different number of elements in each vector, a and b need to have same length')
    end
end

