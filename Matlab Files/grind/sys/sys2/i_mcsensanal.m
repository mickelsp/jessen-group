%this function does an multivariate sensitivity analysis (Klepper 1998) on the  basis of
%the results of a sensitivity analysis. It calculates significance with dummy variables.
%There is no dependence on any GRIND function or variable
%
%mc_output= table with outputs (rows iterations columns model outputs (different variables/times))
%parsets = table with parameter values (rows iterations columns parameters)
%vlabels are labels for the model outputs (separately the output variable and time (nx2 cell matrix));
%plabels are parameter labels
%
%
function [res,a,b]=i_mcsensanal(mc_output, parsets, vlabels, plabels)
if nargout>0
    Soutput = sum_nonan(mc_output);
    Routput = max_nonan(mc_output) - min_nonan(mc_output);
    [res,a,b] = sensmatrix(parsets, mc_output, Soutput,Routput);
    return;
end;
Ndummy = 1000;
Niter = size(mc_output, 1);
Nvars = size(mc_output, 2);
Npars = size(parsets, 2);
Soutput = sum_nonan(mc_output);
Routput = max_nonan(mc_output) - min_nonan(mc_output);
nvlabel =  size(vlabels, 2); %you may have more than one variable label (e.g. var and time)
nplabel =  size(plabels, 2); %you may have more than one parameter label
selectcases = Routput > 1E-12;
if ~all(selectcases)
   warning('GRIND:mcarlo:varsneglected','%d model outputs with negliglible variance removed\n', sum(~selectcases));
   f = find(selectcases);
   mc_output=mc_output(:,f);
   vlabels = vlabels(f, :);
   Nvars = size(mc_output, 2);
   Soutput = sum_nonan(mc_output);
   Routput = max_nonan(mc_output) - min_nonan(mc_output);
end;
%sensitivity matrix
s_matrix = sensmatrix(parsets, mc_output, Soutput,Routput);

%dummy parameters are drawn at random between 0.9-1.1 (not sensitive to this choice).
dummy_matrix = sensmatrix(0.9 + rand(Niter, Ndummy) * 0.2, mc_output, Soutput,Routput);
Lidummy = sqrt(sum(dummy_matrix.^2, 2));
if exist('i_figno','file')
   nfig = i_figno('mcarlo');
else
   nfig = 1;
end;
figure(nfig);
hist(Lidummy, 30);
i_plotdefaults;
title('Distribution total sensitivity of dummy parameters');
xlabel('total sensitivity');
ylabel('frequency');
drawnow;
percentiles = [0.01 0.05 0.1 0.5 0.9, 0.95, 0.99];
signif = i_makepercentiles(Lidummy', percentiles);
percdummy = i_makepercentiles(dummy_matrix', percentiles)';
fprintf('Analysis of %d dummy parameters\n', Ndummy);
for i = 1:length(percentiles)
   if percentiles(i) > 0.5
      fprintf('%g%% percentile: %g\n', percentiles(i) * 100, signif(i));
   end;
end;
% Generate a nice table and copy to clipboard
SensitivityMatrix = cell(Npars + 1+nvlabel+length(percentiles), Nvars + 1+nplabel);
% (2) analyse the sensitivity matrix (statistical toolbox needed?)

Li = sqrt(sum_nonan(s_matrix.^2, 2));
plabels1 = cell(Npars, 1);
plabels2 = plabels1;

for j = 1:Npars
   for i = 1:nplabel
      SensitivityMatrix{j + nvlabel, 1+i} = plabels{j,i};
   end;
   SensitivityMatrix{j +  nvlabel, 1} = Li(j);
   for i = 1:Nvars;
      SensitivityMatrix{j + nvlabel, i + nplabel+1} = s_matrix(j,i);
   end;
   plabels1{j} =  sprintf('%s', plabels{j, :});
   plabels2{j} = sprintf('%5.4f    %s', Li(j), plabels1{j,:});
end;
vlabels1 = cell(Nvars, 1);
for j = 1:nvlabel
   for i = 1:Nvars
      SensitivityMatrix{j, i + nplabel+1} = vlabels{i,j};
      vlabels1{i} = sprintf('%s ', vlabels{i, :});
   end;
end;
for i = 1:length(percentiles)
   SensitivityMatrix{Npars + 1+nvlabel+i, 2} = sprintf('%g%% percentile', percentiles(i) * 100);
   SensitivityMatrix{Npars + 1+nvlabel+i, 1} = signif(i);
end;
for i = 1:Nvars
   for j = 1:length(percentiles)
      SensitivityMatrix{Npars+ 1+nvlabel+j, i + nplabel+1} = percdummy(j,i);
   end;
end;

if exist('uitable','builtin')
   i_table('init', SensitivityMatrix);
else
   varcopy(SensitivityMatrix); %this is the GRIND clipboard function
   disp('The Sensitivity Matrix is copied to the clipboard, paste for instance in Excel')
end;
i_sensplot('init', mc_output,parsets,vlabels1,plabels1);
i_senstimeplot('init', s_matrix,vlabels,plabels1,percdummy);
%plots of each output versus each parameter
%make dendrogram
if ~exist('pdist','file')
   error('GRIND:mcarlo:NoStatToolbox','mcarlo: statistics toolbox required for the cluster analysis');
end;
parndx = i_parseldlg('init', plabels2, Li, percentiles, signif);
%nonans=all(~isnan(scaled));
%scaled = scaled(parndx, nonans);
nonans = all(~isnan(s_matrix));
len=0;
plabels1=plabels1(parndx);
for i=1:length(plabels1)
    if len<length(plabels1{i})
        len=length(plabels1{i});
    end;
end;
plabels2 = plabels2(parndx);
Li=Li(parndx);
for i=1:length(plabels1)
    if len<length(plabels1{i})
        len=length(plabels1{i});
    end;
    plabels2{i}=sprintf('%5.4f %s%s', Li(i),repmat(' ',1,len+1-length(plabels1{i})),plabels1{i});
end;
s_matrix = s_matrix(parndx, nonans);
nlab = length(plabels2);
if nlab<2
    warning('GRIND:mcarlo:dendrogram1var','Cannot make a dendrogram if there is only one parameter/variable');
else
Y = pdist(s_matrix, str2func('sinedist')); % calculate distances
Z = linkage(Y, 'average');
figure(nfig + 3);

dendrogram(Z,nlab,'orientation','right','labels',plabels2);
xlabel('Sine distance');
set(gca, 'fontname','Courier New')
end;

function [s_matrix,regress_a,regress_b] = sensmatrix(p, mc_output, Soutput,Routput)
%create sensitivity matrix
%Sp = sum(p);
%Sp2 = sum(p.^2);
m = sum(~isnan(mc_output),1);
s_matrix = zeros(size(p,2), size(mc_output,2));
if nargout>1
    regress_a=s_matrix;
    regress_b=s_matrix;
end;
pRange = (max(p) - min(p))';
for j = 1:size(p, 2)
   %  pp = repmat(p(:, j), 1, size(mc_output,2)); %parameter values
   % Sxy = sum(pp .* mc_output);
   pp=p(:, j);
   for i = 1: size(mc_output, 2)
         Sxy = sum_nonan(p(:, j) .* mc_output(:, i));
         Sp = sum(pp(~isnan(mc_output(:, i))));
         Sp2 = sum(pp(~isnan(mc_output(:, i))).^2);
         if m(i)>1
            regress = [Sp2 Sp; Sp m(i)] \ [Sxy; Soutput(i)]; %linear regression (conform Grasman)
         else
             regress=[NaN,NaN];
         end;
         a = regress(1);
         if nargout>1
             regress_a(j,i)=regress(1);
             regress_b(j,i)=regress(2);
         end;
         if Routput(i) == 0
            s_matrix(j, i) = NaN; %delen door nul geeft nu NaN (anders zou het extreem gevoelig worden)
         else
            s_matrix(j, i) = a * pRange(j) / Routput(i);
         end;
   end
end

function [d] = sinedist(u, V)
sumX2 = sum_nonan(u.^2, 2);
sumY2 = sum_nonan(V.^2, 2);
sumXY = sum_nonan(V .* repmat(u, size(V, 1), 1), 2);
d = sqrt(1 - (sumXY ./ sqrt(sumX2 .* sumY2)).^2);

function res = sum_nonan(A,p1)
if nargin==1
    p1=1;
end;
A(isnan(A))=0;
res=sum(A,p1);

function res = max_nonan(A,p1,p2)
A(isnan(A))=-inf;
if nargin==1
  res=max(A);
elseif nargin==2
  res=max(A,p1);
elseif nargin==3
  res=max(A,p1,p2);
end;

function res = min_nonan(A,p1,p2)
A(isnan(A))=inf;
if nargin==1
  res=min(A);
elseif nargin==2
  res=min(A,p1);
elseif nargin==3
  res=min(A,p1,p2);
end; 


