function Y = updateSparse_slow(Y,b,indx)
% This is mex-file that updates the values of a sparse matrix Y.
%   The mex-file should be called as:
%
%   updateSparse(Y,b)
%
%   which will implicitly do the following:  Y(omega) = b
%   where "omega" is the set of nonzero indices of Y (in linear
%   ordering, i.e. column-major ordering).
%
%   If "omega" is not sorted, then you must do the following:
%
%       [temp,indx] = sort(omega);  % we don't care about "temp"
%       updateSparse(Y,b,indx);
%
%   which will ensure that everything is in the proper order.

% This file and mex file by Stephen Becker, srbecker@caltech.edu 11.12.08

str = ['Using slow matlab code for updateSparse because mex file not compiled\n',...
    'To disable this warning in the future, run the command:\n',...
    '   warning(''off'',''SVT:NotUsingMex'')\n ',...
    'or install the mex file by running:\n',...
    '   mex updateSparse.c'];

warning('SVT:NotUsingMex',str);

if nargin > 2
    % resort
    b = b(indx);
end

[i,j,s] = find(Y);
[m,n] = size(Y);
Y = sparse(i,j,b,m,n);
