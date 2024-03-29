function [V,EOFs,EC,error]=EOF(D,p)
% function [V,EOFs,EC,error]=EOF2(D,p)
%
% This function decomposes a data set D into its EOFs 
% and eigenvalues V. The most effecient algorithm proposed by
% von Storch and Hannostock (1984) is used.
%
% D = each row is assumed to be a sample. Each column a variable.
% Thus a column represents a time series of one variable
% p = is an optional parameter indicating the number of 
% EOFs the user wants to take into account (smaller than
% the number of time steps and the number of samples).
%
% EOFs= matrix with EOFs in the columns
% V = vector with eigenvalues associated with the EOFs
% EC = EOF Coefficients

% first compute zero-averaged data sets
[n]=size(D,1);
DS = D - repmat(mean(D,1,'omitnan'),n,1);

% Determine size of the data matrix (m=number of time steps, n= number of spatial locations)

q=min(n,p);

% Determine the number of EOFs to be determined (NOE)
if nargin < 2
   NOE=q;
else
   if p=='all'
      NOE=q;
   else
      NOE=min(q,p);
   end
end


% If number of time samples greater than number of variables
% compute the ordinary Covariance Matrix CE
if n>=p    
   CE=DS'*DS/(n-1); % this yields an p x p Covariance Matrix
   [EOFs,S]=eigs(CE,eye(size(CE)),NOE);
end
% If number of time samples smaller than number of variables
% compute the ordinary Covariance Matrix CE
if n < p
   CE = DS*DS'/(n-1); % this yields an n x n Covariance Matrix
   [E,S]=eigs(CE,eye(size(CE)),p);
   EOFs=D'*E;
   for i=1:NOE
      EOFs(:,i)=EOFs(:,i)/norm(EOFs(:,i));
   end
   
end

% rewrite the large eigenvalue matrix to a vector and 
% apply the appropriate normalisation
V=diag(S);

% Determine the EOF coefficients
EC=DS*EOFs;

% Determine the difference between the original data and the 
% reconstructed data
diff=(DS-EC*EOFs');

% determine the L2 error norm for each variable
error=sqrt(sum(abs(diff.^2)));