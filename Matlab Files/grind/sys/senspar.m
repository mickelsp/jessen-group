function senspar(perc, n)
if nargin < 2
   n = 100;
end;
if nargin < 1
   perc = 0.05;
end;
i_senspar('-init',perc,n);
i_senspar('-montecarlo');
i_senspar('-matrix');
i_senspar('-dendrogram');