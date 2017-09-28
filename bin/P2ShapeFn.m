%Authors: Bensingh Dhas and Md Masiur Rahaman
% Description:
%Returns the value of shape functions and its derivatives for Q4 element at any xi and eta.
%Dependencies: None
%Created on: 07th Oct, 2016
function [sp,dsp]=P2ShapeFn(xi)
sp=0.5*[(1-xi)  (xi+1)];
dsp=0.5*[-1 1];
end
