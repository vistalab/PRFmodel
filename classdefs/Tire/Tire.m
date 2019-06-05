classdef Tire < handle
    properties
        car;
        size;
    end
    properties(Dependent, SetAccess = private)
        paintColor;
    end
    
    
    
    
    methods
        function v = Tire
            v.size = [1,2,3];
        end
        function v = get.paintColor(this)
            v = this.car.color;
        end
        
        function plot(this)
            plot(this.size);
        end
    end
end