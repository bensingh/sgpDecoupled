%Authors: Bensingh Dhas and Md Masiur Rahaman
% Description:
%Returns the value of shape functions and its derivatives for Q4 element at any xi and eta.
%Dependencies: None
%Created on: 07th Oct, 2016
function [sp, dsp]=Q4ShapeFn(xi,eta)

sp=0.25*[(xi-1)*(eta-1)  -(xi+1)*(eta-1)  (xi+1)*(eta+1) -(xi-1)*(eta+1)];
dsp=0.25*[(eta-1) -(eta-1) (eta+1) -(eta+1);
    (xi-1) -(xi+1) (xi+1) -(xi-1)];
end
