classdef Car < matlab.mixin.SetGet
    properties(Dependent)
        color;
    end
    properties(GetAccess=public, SetAccess=private)
        tire;
    end
    properties(Access=private)
        bodyColor;
    end
    
    methods
        function set.color(this, v)
            this.bodyColor = v;
        end
        function v = get.color(this)
            v = this.bodyColor;
        end
        function v = getTire(this)
            v = this.tire;
        end
        function obj = Car()
            obj.tire = Tire;
            obj.tire.car = obj;
            obj.color = 'black';
        end
    end
end
