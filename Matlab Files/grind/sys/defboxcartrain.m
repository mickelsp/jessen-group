function boxnr = defboxcartrain(name)
global g_grind;
if isfield(g_grind, 'boxcar')
   boxnr = 1;
   while boxnr <= length(g_grind.boxcar.names) && ~strcmp(g_grind.boxcar.names{boxnr}, name)
      boxnr = boxnr + 1;
   end;
   if boxnr > length(g_grind.boxcar.names)
      %      boxnr=[];
      %   end;
      %   boxnr = strmatch(name, g_grind.boxcar.names, 'exact');
      %  if isempty(boxnr)
      boxnr = length(g_grind.boxcar.trains) + 1;
      g_grind.boxcar.names = [g_grind.boxcar.names {name}];
   end;
else
   boxnr = 1;
   g_grind.boxcar.names{1} = name;
end;
boxcartrain.name = name;
boxcartrain.gcycl = [];
boxcartrain.outflow = 0;
boxcartrain.flow = [];
g_grind.boxcar.trains{boxnr} = boxcartrain;
