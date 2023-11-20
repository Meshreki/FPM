function [ O,P ] = GDUpdate_Multiplication_rank1(O,P,dpsi,Omax,cen,Ps,alpha,beta,step_size)
%GDUPDATE_MULTIPLICATION update estimate of O and P according to gradient
%descent method, where psi = O*P
%   Inputs:
%   O0: object estimate, n1xn2
%   P0: pupil function estimate: m1xm2
%   psi: update estimate field estimate
%   psi0: previous field estimate
%   cen: location of pupil function
%   alpha: gradient descent step size for O
%   betta: gradient descent step size for P
%   Ps: support constraint for P0, e.g. spatially confined probe or
%   objective with known NA
%   iters: # of iterations to run on updates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All main.m, AtlerMin.m, GDUpdate_Multiplication_rank1.m, USAF_Parameter.m, 
% hela_Parameter.m files are modified and extended to:
%
% 1- allowing a FULL reconstruction of the input images for 
% any dimesnions (i.e. not only square crops of them!).
% 2- saving images (and figures) for multiple variables.
% 3- saving dirac peaks positions.
% 4- adding more descriptive comments (in addition to those added by Ivo Ihrke).
% 4- investigating the huge error of the Algorithm to be due to pixels w/ 
% very large (intensity) values in the input stacks (and not due to
% the Algorithm itself!).
%
% last modified on 27.05.2022
% by John Meshreki, john.meshreki@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% last modified by Lei Tian, lei_tian@alum.mit.edu, 3/1/2014


% size of P, Np<=No
Np = size(P); Np = Np(:); % convert a row into a column

% operator to put P at proper location at the O plane
n1 = [cen(1)-floor(Np(1)/2), cen(2)-floor(Np(2)/2)];
n2 = [n1(1)+Np(1)-1, n1(2)+Np(2)-1];
% operator to crop region of O from proper location at the O plane
%downsamp = @(x) x(n1(1):n2(1),n1(2):n2(2));
downsamp = @(x) x(n1(1)+1:n2(1)+1,n1(2)+1:n2(2)+1);

O1 = downsamp(O);

%Fig. 3 O-and P-updates for 1 LED
O(n1(1):n2(1),n1(2):n2(2)) = O(n1(1):n2(1),n1(2):n2(2))...
    + step_size * 1/max(max(abs(P)))*abs(P).*conj(P).*dpsi./(abs(P).^2+alpha);
P = P+1/Omax*(abs(O1).*conj(O1)).*dpsi./(abs(O1).^2+beta).*Ps;

end

