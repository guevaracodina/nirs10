function M = Mausphercnst(X,alpha,verbose)
%MAUSPHERCNST Mauchly's sphericity test for orthogonal contrasts.
% Other hypothesis used with the Mauchly's sphericity test is on the orthonormalized
% contrasts
%
% H0: C'SC = s^2I  ,where C is any full-rank (p-1) x p matrix of orthonormal
% contrasts.
%
% A is the nonorthonormalized contrasts matrix
%
% C = GMO(A') is a (p-1) × p column-orthonormal matrix (by the Gram-Schmidt 
% orthonormalization procedure)
%
% L = q^q|C'SC|/[tr(C'SC)]^q  
%
% LL = - (N - 1 - ((2q^2 + q + 2)/(6q)) ln L
%
%    = - (N - 1 - ((2q^2 + q + 2)/(6q)) [ln|C'SC|- q ln [tr(C'SC)/q]]

% This statistic has an approximate Chi-square distribution with 0.5q(q + 1) - 1 =
% q(q + 1)/2 - 1 degrees of freedom.
%
%     Syntax: function [Mausphercnst] = Mausphercnst(X,alpha) 
%      
%     Inputs:
%          X - multivariate data matrix. 
%      alpha - significance level (default = 0.05). 
%     Output:
%          n - sample-size
%          p - variables
%          L - Mauchly's statistic used to test any deviation from
%              an expected sphericity
%         df - degrees of freedom
%         X2 - Chi-square statistic
%          P - probability that null Ho: is true
%
% Example: From the example in Diggle, Liang and Zeger (1994, Table 3.2, p. 35). Body weights
%        of 48 pigs measured in 9 consecutive weeks of followup. With a significance level of
%        0.05, we are interested to test the sphericity taken into account the constrasts
%        approach.
%
% Total data matrix must be:
%X=[24.0 32.0 39.0 42.5 48.0 54.5 61.0 65.0 72.0;22.5 30.5 40.5 45.0 51.0 58.5 64.0 72.0 78.0;
%   22.5 28.0 36.5 41.0 47.5 55.0 61.0 68.0 76.0;24.0 31.5 39.5 44.5 51.0 56.0 59.5 64.0 67.0;
%   24.5 31.5 37.0 42.5 48.0 54.0 58.0 63.0 65.5;23.0 30.0 35.5 41.0 48.0 51.5 56.5 63.5 69.5;
%   22.5 28.5 36.0 43.5 47.0 53.5 59.5 67.5 73.5;23.5 30.5 38.0 41.0 48.5 55.0 59.5 66.5 73.0;
%   20.0 27.5 33.0 39.0 43.5 49.0 54.5 59.5 65.0;25.5 32.5 39.5 47.0 53.0 58.5 63.0 69.5 76.0;
%   24.5 31.0 40.5 46.0 51.5 57.0 62.5 69.5 76.0;24.0 29.0 39.0 44.0 50.5 57.0 61.5 68.0 73.5;
%   23.5 30.5 36.5 42.0 47.0 55.0 59.0 65.5 73.0;21.5 30.5 37.0 42.5 48.0 52.5 58.5 63.0 69.5;
%   25.0 32.0 38.5 44.0 51.0 59.0 66.0 75.5 86.0;21.5 28.5 34.0 39.5 45.0 51.0 58.0 64.5 72.5;
%   31.0 38.0 48.0 54.0 60.0 62.0 66.5 75.5 84.0;27.5 32.5 36.0 43.0 49.5 52.5 56.0 61.0 64.0;
%   30.0 37.0 45.0 51.0 58.0 63.0 67.5 74.5 81.0;26.0 32.0 40.5 45.5 52.5 55.5 62.5 69.5 74.0;
%   26.0 32.5 39.5 44.0 48.0 54.5 58.0 66.0 73.0;28.5 35.5 41.5 47.5 54.0 59.5 63.5 71.0 78.5;
%   26.5 34.5 42.0 48.5 55.5 62.0 68.0 76.5 85.0;27.5 33.5 41.0 45.0 50.5 56.0 62.5 71.0 78.0;
%   22.5 27.0 33.5 38.5 41.0 49.0 56.0 64.0 68.0;22.0 26.5 32.5 38.5 43.5 50.5 56.5 63.5 68.5;
%   23.5 29.0 35.5 40.0 45.0 50.0 56.5 63.0 67.5;22.5 29.5 36.5 42.0 45.0 55.0 61.0 68.0 72.0;
%   27.5 34.5 42.0 47.5 53.0 63.0 72.0 79.0 85.5;23.5 28.0 33.0 37.0 38.5 48.0 52.5 62.0 64.5;
%   24.5 30.0 38.5 42.0 47.5 54.0 62.5 71.5 77.0;24.5 31.5 40.5 46.5 51.5 61.5 68.5 77.5 84.5;
%   24.5 32.0 39.0 45.0 51.0 55.5 61.5 69.0 75.5;24.0 32.5 40.0 48.0 54.5 61.5 68.0 74.5 81.0;
%   24.0 31.5 38.5 44.0 51.5 57.5 64.0 72.5 79.0;24.5 32.5 39.5 44.5 52.5 56.5 62.0 67.5 72.5;
%   24.5 32.0 38.5 44.0 50.0 56.0 63.5 69.5 76.0;25.5 33.0 41.5 47.0 55.5 60.5 66.5 77.0 82.0;
%   25.5 32.0 39.0 45.5 51.0 57.5 63.5 72.0 78.5;25.0 31.0 36.5 43.0 50.5 55.0 62.5 69.0 75.5;
%   26.5 30.5 33.0 39.0 43.5 49.5 56.5 61.0 65.0;24.0 32.0 39.0 44.5 50.0 56.0 63.0 67.5 74.0;
%   24.5 31.0 37.5 43.5 48.0 56.0 62.5 66.5 70.5;27.0 34.5 42.0 48.5 53.0 60.0 67.0 73.0 76.0;
%   31.0 39.0 47.5 51.0 57.0 64.0 71.0 77.0 80.5;27.0 33.5 40.0 46.5 53.0 60.0 66.5 72.5 80.0;
%   29.5 37.0 46.0 52.5 60.0 67.5 76.0 81.5 88.0;28.5 36.0 42.5 49.0 55.0 63.5 72.0 78.5 85.5];
%
% Calling on Matlab the function: 
%     Mausphercnst(X)
%
%             Answer is:
%  --------------------------------------------------------------------------
%   Sample-size    Variables        L         df        Chi-square       P
%  --------------------------------------------------------------------------
%        48            9          0.0001      35         430.5400      0.0000
%  --------------------------------------------------------------------------
%  With a given significance level of: 0.05
%  Assumption of sphericity is not tenable.
%
%  Created by A. Trujillo-Ortiz, R. Hernandez-Walls and K. Barba-Rojo
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.mx
%
%  Copyright. April 18, 2008.
%
% To cite this file, this would be an appropriate format:
% Trujillo-Ortiz, A., R. Hernandez-Walls and K. Barba-Rojo. (2008). Mausphercnst:
%   Mauchly's sphericity test for orthogonal contrasts. A MATLAB file. [WWW document].
%   URL http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=19648
%
%  References:
%
%  Davies, A. W. (1971), Percentile approximations for a class of likelihood ratio
%              criteria. Biometrika, 58:349-356.
%  Diggle, P. J., Liang, K. and and Zeger, S. L. Zeger. (1994), Analysis of 
%              Longitudinal Data. Oxford Science Publications. 
%  Johnson, R. A. and Wichern, D. W. (1992), Applied Multivariate Statistical Analysis.
%              3rd. ed. New-Jersey:Prentice Hall. pp. 158-160.
%  Mauchly, J. W. (1940), Significance test for sphericity of a normal n-variate
%              distribution. The Annals of Mathematical Statistics, 11:204-209.
%
  
if nargin < 2, 
   alpha = 0.05;  %(default)
end; 

if nargin < 1, 
   error('Requires at least one input arguments.');
end;

[n,p] = size(X);
S = cov(X);  %variance-covariance matrix
R = diag(ones(p,1))-diag(ones(p-1,1),1);
q = p-1;
A = R(1:q,:);  %nonorthonormalized contrasts matrix
B = A';
C = cgrscho(B);  %m×n column-orthonormal output matrix
L = q^q*det(C'*S*C)/(trace(C'*S*C)^q);  %Mauchly's statistic
LL = -(n-1-((2*q^2+q+2)/(6*q)))*log(L);  %approximate to Chi-square distribution
F = q*(q+1)/2-1;  %degrees of freedom
P = 1-chi2cdf(LL,F);  %Chi-square probability-value

M.n = n;
M.p = p;
M.L = L;
M.df = F;
M.X2 = LL;
M.P = P;
if verbose
fprintf('--------------------------------------------------------------------------\n');
disp(' Sample-size    Variables        L         df        Chi-square       P')
fprintf('--------------------------------------------------------------------------\n');
fprintf('%8.i%13.i%16.4f%8.i%17.4f%12.4f\n',n,p,L,F,LL,P);
fprintf('--------------------------------------------------------------------------\n');
fprintf('With a given significance level of: %.2f\n', alpha);
end     
if P >= alpha;
   M.H = 0; 
   if verbose
   fprintf('Assumption of sphericity is tenable.\n\n');
   end
else
   M.H = 1;
   if verbose
   fprintf('Assumption of sphericity is not tenable.\n\n');
   end
end 

return;