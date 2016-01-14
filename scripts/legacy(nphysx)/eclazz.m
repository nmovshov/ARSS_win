function eqcls=eclazz(equiv,data)
%ECLAZZ   Partition a set into equivalence classes using a given relation
% ECLAZZ is an m-language version of the Numerical Recipes function eclazz.
% The equivalence relation equiv divides some dataset data into
% equivalence classes (subsets). Recall that an equivalence relation is any
% relation that satisfies reflexivity, symetry, and transitivity.
% Intuitively this is a general "same-as" relation. The task is to use this
% relation to efficiently construct the equivalence subsets of a given set.
%
% Numerical Recipes (3rd edition) contains an efficient non-recursive
% algorithm due to D. Eardley [1] in c++. This m-function follows the same
% algorithm.
%
% Inputs:
%  equiv: function handle to a valid RST function of two integer
%   indices. equiv(i,j) will test the equivalence of the i-th and j-th
%   elements of data. **No checks will be made on the validity of equiv.**
%   data: 1D cell array containing the set to be partitioned.
% Oututs:
%  eqcls: 1D array of integer labels. eqcls(i) is the integer label of the
%  equivalence class that contains data{i}.
%
% Author: Naor Movshovitz.
% References: [1] Numerical Recipes (3rd edition) ch 8.6.

%% minimal input filter
error(nargchk(2,2,nargin))
if ~isa(equiv,'function_handle'), error('don''t be a dick'), end
if ~isscalar(equiv), error('don''t be a dick'), end
if ~iscell(data), error('don''t be a dick'), end
if ~isvector(data), error('don''t be a dick'), end

%% initialize output array for speed
%eqcls=uint32(zeros(size(data))); % int classes don't support NaN
eqcls=nan(size(data)); % is uint32 better or worse or nothing?

%% The NR algorithm (with the book's exact, and useless, comments)
eqcls(1)=1;
for j=2:length(data) % Loop over first element of all pairs.
    eqcls(j)=j;
    for k=1:j-1      % Loop over second element of all pairs.
        eqcls(k)=eqcls(eqcls(k)); % Sweep it up this much.
        if equiv(data{j},data{k}), eqcls(eqcls(eqcls(k)))=j; end % Good exercise for the reader to figure out why this much ancestry is necessary!
    end
end
for j=1:length(data)
    eqcls(j)=eqcls(eqcls(j)); % Only this much sweeping is needed finally.
end