function [varlist] = i_savedata(varlist, amatrix, checkonly)
global g_grind g_data;
if nargin < 3
   checkonly = 0;
end;
isdata = 0;
hast = 0;
%err = 0;
for i = 1:size(varlist, 2)
   if ~isempty(varlist{i})&&~strcmp(varlist{i},'###skip###');
      k = i_getno(varlist{i});
      if isempty(k.no) && (strcmp(varlist{i}, 't'))
         if hast
 %           err = 1;
            error('GRIND:setdata:TwoTimeCols','two columns represent time (t)');
         end;
         hast = 1;
         for i1=2:size(amatrix,1)
            if amatrix(i1)<amatrix(i1-1)
 %              err=1;
               error('GRIND:setdata:DecrTime','Error in setdata: time must be monotonously increasing');
            end;
         end;
       elseif ~isempty(k.no) && k.ispar
         s = sprintf('Column %d (%s) represents a parameter, do you want to make it an external variable?' ...
            , i, varlist{i});
         button=questdlg(s,'data','Yes','Skip','Rename','Yes');
         switch button
          case 'Yes'
            answer=inputdlg({sprintf('What is the default value of %s?',varlist{i})},'data',1,{'0'});
            defextern(varlist{i}, answer);
          case 'Skip'
            varlist{i} = [];
          case 'Rename'
            answer=inputdlg(sprintf('What is the new name of column %d?',i),'data',1,varlist(i));
            if ~isempty(answer)
               varlist{i} = answer;
            end;
         end;
      elseif isempty(k.no)
         l = findextern(varlist{i});
         if isempty(l)
            s = sprintf('The name of column %d (%s) is unknown, do you want to rename it?' ...
               ,i, varlist{i});
            button=questdlg(s,'data','Yes','Skip','Yes'); 
            if strcmp(button, 'Yes')
               answer=inputdlg({sprintf('What is the new name of column %d?',i)},'data',1,varlist(i));
               if ~isempty(answer)
                  varlist{i} = answer{1};
               end;
            else
               varlist{i} = [];
            end;
         end;
      end;
   end;
end;
for i = 1:size(varlist, 2)
   if ~isempty(varlist{i})
      k = i_getno(varlist{i});
      if ~isempty(k.no) && ~k.ispar && ~k.isfun
         isdata = 1;
      end;
   end;
end;
if ~checkonly
   if isdata
      if ~isfield(g_data, 'options')
         opts = optimset('fminsearch');
         g_data.options.Display = opts.Display;
         g_data.options.TolX = 1E-6;
         g_data.options.TolFun = 1E-6;
         g_data.options.MaxFunEvals1 = '200*numberOfVariables';
         g_data.options.MaxIter1 = '200*numberOfVariables';
         g_data.options.PosOnly = 0;
         g_data.options.ResetEachStep = 0;
      end;
      g_data.obs = NaN+zeros(size(amatrix, 1), g_grind.statevars.dim);
      g_data.pred = zeros(size(amatrix, 1), g_grind.statevars.dim);
   end;
   datat = [];
   for i = 1:size(varlist, 2)
      k = i_varno(varlist{i});
      if isempty(k) && (strcmp(varlist{i}, 't'))
         datat = amatrix(:, i);
      elseif isempty(k)
         l = findextern(varlist{i});
         if ~isempty(l)
            if isempty(datat)
               g_grind.externvars{l}.data = amatrix(:, i);
            else
               g_grind.externvars{l}.data = [datat, amatrix(:,i)];
            end;
            g_grind.externvars{l}.data=g_grind.externvars{l}.data(~isnan(sum(g_grind.externvars{l}.data,2)),:);
            g_grind.lastsettings={};
         end;
      elseif ~isempty(k)
         g_data.obs(:, k) = amatrix(:, i);
      end;
   end;
   if isdata
      g_data.t = datat;
      g_data.minobs = min(g_data.obs);
      g_data.maxobs = max(g_data.obs);
   end;
end;
function varno = findextern(s)
global g_grind;
varno = [];
for k = 1:length(g_grind.externvars)
   if strcmp(s, char(g_grind.externvars{k}.name))
      varno = k;
      return;
   end;
end;

