classdef zernike
% -------------------------------------------------------------------------
% Copyright C 2016 Amir Tahmasbi
% Texas A&M University
% amir.tahmasbi@tamu.edu
% http://www.amirtahmasbi.com/

% License Agreement: To acknowledge the use of the code please cite the 
%                    following papers:
%
% [1] A. Tahmasbi, F. Saki, S. B. Shokouhi, 
%     Classification of Benign and Malignant Masses Based on Zernike Moments, 
%     Comput. Biol. Med., vol. 41, no. 8, pp. 726-735, 2011.
%
% [2] F. Saki, A. Tahmasbi, H. Soltanian-Zadeh, S. B. Shokouhi,
%     Fast opposite weight learning rules with application in breast cancer 
%     diagnosis, Comput. Biol. Med., vol. 43, no. 1, pp. 32-41, 2013.
% 
% Purpose:
%   Claculate Zernike polynomials and moments of order n and repetition m 
%
% Class properties:
%   order       (non-negative scalar)
%   repetition  (integer such that order - abs(repetition) = even)
%   rho         (scalar radial distance 0<= rho <=1)
%
% Revision history:
%   7/13/2016 - Created from existing code.
%
% -------------------------------------------------------------------------

   properties (Access = public)
       order;               % the order of the Zernike polynomial
       repetition;          % the repetition of the Zernike polynomial
       rho;                 % radial distance 0<= rho <=1
   end
   
   methods
       % class constructor
        function obj = zernike(varargin)
           if nargin>0
               if (~mod(varargin{1},1) && varargin{1}>=0)
                    obj.order = varargin{1};
               else
                    error('The order must be a non-negative integer.');
               end
               if ~mod(varargin{1} - abs(varargin{2}), 2)
                    obj.repetition = varargin{2};
               else
                    error('The difference between the order and the repetition must be an even number.');
               end
           end 
        end % end constructor 
       
       % get radial polynomial
        function [output, obj] = getRadialPoly(obj, rho)
            % check for any errors
            if nargin>0
                if isvector(rho)
                    obj.rho = rho(:);
                else
                    error('The radial distance must be a vector.');
                end
            end
            % obtain the order and repetition
            n = obj.order;
            m = obj.repetition;
            output = zeros(size(obj.rho));      % initilization
            
            % compute the radial polynomial
            for s = 0:(n-abs(m))/2
              c = (-1)^s*factorial(n-s) / ...
                (factorial(s)*factorial((n+abs(m))/2-s)*factorial((n-abs(m))/2-s));
              output = output + c * obj.rho.^(n-2*s);
            end  
        end % end getRadialPoly method
       
   end
end