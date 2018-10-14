function make

% Recompile all mex functions in the tools directory
%
% Christian Merkwirth 2002

clear functions

mymex('nn_search.cpp', '-DBRUTE', '-output', 'nn_brute')
mymex('nn_search.cpp', '-DPARTIAL_SEARCH', '-output', 'nn_search')
mymex('nn_prepare.cpp', '-DPARTIAL_SEARCH', '-output', 'nn_prepare')
mymex('emb_nn_search.cpp', '-DPARTIAL_SEARCH' ,'-DBRUTE', '-DMAX_NORM','-output', 'emb_brute_nn')
mymex('range_search.cpp', '-DBRUTE', '-output', 'range_brute')
mymex('range_search.cpp', '-DPARTIAL_SEARCH')
mymex('corrsum2.cpp', '-DPARTIAL_SEARCH')
mymex('corrsum.cpp', '-DPARTIAL_SEARCH')

function mymex(target, varargin)

if isunix
	mex('-v', '-I..', '-I.', varargin{:},  '-O', target)
else
	mex('-v', '-I..', '-I.', varargin{:}, '-O', target) 
end
