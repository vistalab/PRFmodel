function cl_out = pmFindAttributes(obj,attrName,varargin)
   if ischar(obj)
      mc = meta.class.fromName(obj);
   elseif isobject(obj)
      mc = metaclass(obj);
   end
   ii = 0; numb_props = length(mc.PropertyList);
   cl_array = cell(1,numb_props);
   for  c = 1:numb_props
      mp = mc.PropertyList(c);
      if isempty (findprop(mp,attrName))
         error('Not a valid attribute name')
      end
      attrValue = mp.(attrName);
      if attrValue
         if islogical(attrValue) || strcmp(varargin{1},attrValue)
            ii = ii + 1;
            cl_array(ii) = {mp.Name};
         end
      end
   end
   cl_out = cl_array(1:ii);
end