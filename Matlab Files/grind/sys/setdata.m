%SETDATA   Enter/edit data for external or state variables
%   Command to enter data for parameter optimizing or for external variables. 
%   You can either enter the data manually or as a file. (see for formats: 
%   loaddata)
%
%   Usage:
%   SETDATA - Enter/edit data matrix
%   SETDATA AMatrix - Use AMatrix as data matrix.
%   SETDATA(aMatrix,varlist) aMatrix=matrix with data
%   varlist= list with variables (t = time) (if omitted,
%   it is assumed that the first column is time and the
%   other columns are in same order as g_grind.statevars).
%   SETDATA -CLEAR - clear the current matrix
%  
%    See also loaddata, optimpars, defextern, g_data 

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function setdata(amatrix, varlist)%,TolX,TolFun,MaxFunEvals,MaxIter)
global g_data g_grind;
i_parcheck;
if nargin == 1
   if ischar(amatrix) && strcmpi(amatrix, '-clear')
      g_data.obs = [];
      g_data.pred = [];
      g_data.t = [];
      g_data.pars = [];
      g_data.minobs = [];
      g_data.maxobs = [];
      g_grind.lastsettings = [];
      out '-silent' '-remove' '-data';
      out '-silent' '-remove' '-extern';
      for i = 1:length(g_grind.externvars)
         g_grind.externvars{i}.data =  [];
      end;
      return;
   end;
end;
if nargin  > 0
   if ischar(amatrix)
      amatrix = i_checkstr(amatrix);
   elseif iscell(amatrix)
      amatrix = i_checkstr(char(amatrix));
   end;
end;
TAB = sprintf('\t');
if nargin == 1
   disp('Assumed order of columns:');
   varlist=cell(1,length(g_grind.externvars));
   for i = 1:length(g_grind.externvars)
      varlist{i} = g_grind.externvars{i}.name;
   end;
   varlist = [{'t'} g_grind.statevars.names varlist];
elseif nargin == 0
   amatrix = [];
   tr = [];
   if ~isempty(g_data);
      if isfield(g_data, 'obs')
         amatrix = g_data.obs;
      end;
      if isfield(g_data, 't')
         tr = g_data.t;
      else
         tr = [];
      end;
   end;
   for i = 1:length(g_grind.externvars)
      if ~isempty(g_grind.externvars{i}.data)
         [tr, amatrix] = i_concatdata(tr, amatrix, g_grind.externvars{i}.data(:, 1), g_grind.externvars{i}.data(:, 2));
      end;
   end;
   varlist=cell(1,length(g_grind.externvars));
   for i = 1:length(g_grind.externvars)
      varlist{i} = g_grind.externvars{i}.name;
   end;
   varlist = [{'t'} g_grind.statevars.names varlist];

   [amatrix,varlist]=i_setdatadlg('init',[tr,amatrix],varlist);
elseif nargin > 1
   if ischar(varlist)
      l = length(varlist) + 1;
      l2 = l - 1;
      while l2 < l
         varlist = strrep(varlist, [TAB TAB], [TAB '###skip###' TAB]);
         varlist = strrep(varlist, TAB, ' ');
         varlist=strrep(varlist,'  ',' ');
         l = l2;
         l2 = length(varlist);
      end;
      varlist=['{''' strrep(strtrim(varlist),' ',''' ''') '''}'];
      varlist = memo2str(varlist, 0);
   end;
   varlist = i_checkstr(varlist);
end;
i_savedata(varlist, amatrix);
out('-silent','-add','-data','-extern');


function s = memo2str(answer, addmatrixthings)
s1 = i_memo2cell(answer);
%s = '';
f=strfind(s1{1}, '=');
if ~isempty(f)
   s1{1} = s1{1}(f(1) + 1:length(s1{1}));
end;
for i = 1:length(s1)
   f = strfind(s1{i},'...');
   if ~isempty(f)
      s1{i} = s1{i}(1:f(1) - 1);
   end;
    if addmatrixthings
      f = strfind(s1{i}, ';');
       if isempty(f)
          s1{i}=[s1{i} ';'];
       end;
    end;   
%    if addmatrixthings
%       f = strfind(s1{i}, ';');
%       if isempty(f)
%          s = [s ';' s1{i}];
%       else
%          s = [s s1{i}];
%       end;
%    else
%       s = [s s1{i}];
%    end;
end;
s=sprintf('%s',s1{:});
if addmatrixthings && isempty(strfind(s, '['))
   s=['[' s ']'];
end;


