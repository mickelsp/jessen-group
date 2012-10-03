% CREATE_NETWORK - create different kinds of networks
% create a matrix with ones for connections and zeros for not, using different rules:
%    -exponential (e)
%    -scale-free (s)
%    -lattice 4 neighbors (l)
%    -random (r) pure (!) random, no check for unconnected nodes
%
%  Usage:
%    [A,ndx]=CREATE_NETWORK(N, TYPE, PAR) N= number of nodes, type is kind of network ('e','s','l' etc)
%       PAR = a specific parameter (see below). The function returns A (sparse matrix with connections (use FULL to convert
%       to full matrix) and ndx, which  are the indices>0 of the matrix A
%
%    PAR: for exponential and scale-free PAR = the number of legs
%         for lattice PAR = 1 for periodic boundaries (and 0 for non periodic boundaries)
%         for random PAR = probability of a connection
%
%  See also:
%     SPARSE2NDX, NDX2SPARSE
%
function [A, ndxs] = create_network(n, type, par1)

if nargin < 1
   prompt={'Number of nodes','Kind of network (exponential (e), scale-free (s), lattice (l) or random (r)','Extra parameter (number of legs/probability)'};
   def={'100','e','1'};
   dlgTitle = 'Create Network';
   lineNo = 1;
   answer = inputdlg(prompt, dlgTitle, lineNo, def);
   n = str2double(answer{1});
   type = answer{2};
   par1 = str2double(answer{3});
elseif nargin < 2
   type = 'e';
   disp('exponential network'); %default network
end;
random_order = 0; % switch to add the nodes in random order or not
if ~ischar(type)
   error('GRIND:create_network:UnknownType','Network type "%g" not supported', type);
end;

if strncmpi(type, 'e', 1)
   %%%% EXPONENTIAL NETWORK %%%%
   if nargin < 3
      legs = 1; %default number of legs is one
   else
      legs = round(par1); %par1 are the number of legs
      if legs < 1
         error('GRIND:create_network:TooFewLegs','The number of new connections ("legs") should be at least 1');
      end;
   end;
   if random_order
      %shuffle the indexes of nodes (1:n in random order)
      rands = rand(n, 1); %#ok
      [x, ndx] = sort(rands); %ndx are the nodes in random order
   else
      ndx = (1:n)';
   end;
   ndx1 = ndx;
   %the nodes in random order are added one by one
   %now set for each node the number of nodes that should be already in the network
   maxndx = (0:n - 1)';
   ndxs = [];
   %draw the nodes to connect
   for j = 1:legs
      %choose a random node that is already in the network (don't mind double connections)
      ndx2 = ndx1(floor(rand(n, 1) .* maxndx) + 1);
      %save the subscripts of the connections
      ndxs = [ndxs; [ndx1(2:n), ndx2(2:n)]]; %no check for doubles
   end;
   
   %fill the square matrix with connections
   A=ndx2sparse(ndxs);
   
   
elseif strncmpi(type, 's', 1)
   
   %%%% SCALE-FREE NETWORK %%%%
   % is similar as exponential, but with weighted selection
   if nargin < 3
      legs = 1; %default number of legs is one
   else
      legs = round(par1); %par1 are the number of legs
      if legs < 1
         error('GRIND:create_network:TooFewLegs','The number of new connections ("legs") should be at least 1');
      end;
   end;
   if random_order
      %shuffle the indexes of nodes (1:n in random order)
      rands = rand(n, 1); %#ok
      [x, ndx] = sort(rands); %ndx are the nodes in random order
   else
      ndx = (1:n)';
   end;
   %the nodes in random order are added one by one
   %now set for each node the number of nodes that should be already in the network
   counts = zeros(n, 1);  % count the links
%   A = sparse(0);
   ndxs = zeros(n * legs, 2);
   k = 1;
   %draw the nodes to connect
   for i = 2:n
      if (i == 2)  %first node is added
         ndxs(k, :) = ndx([i, i - 1])'; %ndx could be shuffled
         k = k + 1;
         counts([i, i - 1]) = counts([i, i - 1]) + 1;
      else
         P = counts ./ sum(counts);
         P = P(1:i - 1);
         for j = 1:legs % could be made more efficient if length(added) > legs
            added = [];
            while isempty(added)
               r = rand(i - 1, 1);
               added = find(r < P);
            end;
            anndx = [i, added(floor(rand(1) .* length(added)) + 1)]; %if more than one added,
            %choose at random one of these
            %(THIS IS DIFFERENT FROM LUISJO, he chooses allways the first);
            counts(anndx) = counts(anndx) + 1;
            ndxs(k, :) = ndx(anndx)';
            k = k + 1;
         end;
      end;
   end;
   %fill the square matrix with connections
   A=ndx2sparse(ndxs);
   
elseif strncmpi(type, 'r', 1)
   
   %%%%% RANDOM NETWORK %%%%%
   %pure random no check if nodes are all connected
   if nargin < 3
      prob = 0.5; %default probability of a connection = 0.5
   else
      prob = par1; %par1 is the probability
   end;
   A = triu(rand(n)); %we only need one triangle as the matrix should be symmetric
   A(logical(speye(size(A)))) = 0;
   A = sparse(A > 1 - prob);  %we have a connection with a probability of prob (1 - p > rand)
   A =A+triu(A)'; %make symmetric matrix
   if nargout > 1
      ndxs = sparse2ndx(A);
   end;
   
elseif strncmpi(type, 'l', 1)
   %%%%% REGULAR LATTICE (4 neighbor) %%%%%
   % n can now be a vector [sizeX,sizeY]
   % if not n is assumed to be the total number of cells of a square matrix of [sqrt(n),sqrt(n)]
   if length(n) == 1
      siz = [round(sqrt(n)), round(sqrt(n))];
      if sqrt(n) - siz(1) > 0.01
         warning('GRIND:create_network:square','Matrix size assumed to be square: changed to [%d,%d] ',siz);
      end;
   else
      siz = n;
   end
   if nargin < 3
      isround = 1; %default periodic boundaries
   else
      isround = par1; %par1 is the probability
   end;
   ntot=prod(siz);
   ndxs = zeros(4*ntot,2);
   k=1;
   for i = 1:siz(1)
      for j = 1:siz(2);
         ii = sub2ind(siz, i, j);
         ndxs(k:k+3,1)=ii;
         ndxs(k,2)= sub2ind(siz, i, getndx(j - 1, siz(2), isround));
         ndxs(k+1,2)=sub2ind(siz, i, getndx(j + 1, siz(2), isround));
         ndxs(k+2,2)=sub2ind(siz, getndx(i + 1, siz(1), isround), j);
         ndxs(k+3,2)=sub2ind(siz, getndx(i - 1, siz(1), isround), j);
         %neighbors getndx takes the periodic boundaries into account
         k=k+4;
       end;
   end;
   A=ndx2sparse(ndxs);
  % A=ndx2sparse(ndxs(ndxs(:,1)>0,:));
else
   error('GRIND:create_network:UnknownType','Network type "%s" not supported', type);
end;

function i = getndx(i, s1, isround)
if isround
   if i <= 0
      i = i + s1;
   elseif i > s1
      i = i - s1;
   end;
else
   if i <= 0
      i = 1;
   elseif i > s1
      i = s1;
   end;
end;

