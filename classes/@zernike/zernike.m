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
%   Calculate Zernike polynomials and moments of order n and repetition m 
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

%TODO Zernike transform and its inverse in a separate class

   properties (Access = public)
       order;               % the order of the Zernike polynomial
       repetition;          % the repetition of the Zernike polynomial
       rho;                 % radial distance 0<= rho <=1
       noll_index;          % Noll index of the moment
       aberration_name;     % classical name of the aberration
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
                    obj = setAberrationTypeUtil(obj);
               else
                    error('The difference between the order and the repetition must be an even number.');
               end
           end 
        end % end constructor 
       
        % -------------------------
        % get the radial polynomial
        function [output, obj] = getRadialPoly(obj, rho)
            % check for any errors
            if nargin>1
                obj.rho = rho;
            end
            % obtain the order and repetition
            n = obj.order;
            m = obj.repetition;
            output = zeros(size(obj.rho));      % initilization
            
            % compute the radial polynomial
            for s = 0:(n-abs(m))/2
              c = (-1)^s*factorial(n-s) / ...
                (factorial(s)*factorial((n+abs(m))/2-s)*factorial((n-abs(m))/2-s));
              output = output + c * obj.rho .^ (n-2*s);
            end  
        end % end getRadialPoly method
        
        % ----------------------
        % get the Zernike moment
        function [output, pupil, obj] = getZernikeMoment(obj, img)
            
            [N, M]  = size(img);   
            x       = -1+1/M:2/M:1-1/M; 
            y       = -1+1/N:2/N:1-1/N;
            [X,Y]   = meshgrid(x,y);
            
            [theta, r]  = cart2pol(X, Y);         
            obj.rho     = (r<=1) .* r;
            
            R = obj.getRadialPoly();        % get the radial polynomial
            
            if obj.repetition >= 0          
                pupil = R .* cos(obj.repetition * theta);
            else
                pupil = R .* sin(-obj.repetition * theta);
            end

            Product = img .* pupil;
            output = sum(Product(:));       % calculate the moments

            cnt = nnz(R)+1;                 % count the number of pixels inside the unit circle
            output = (obj.order+1)*output/cnt; % normalize the amplitude of moments
            
        end % end getZernikeMoments method
        
        % -----------------------
        % get the aberration type
        function output = getAberrationType(obj)
            output = obj.aberration_name;
        end
        
        % -----------------------
        % set the aberration type 
        function obj = setAberrationType(obj, type)
            switch type
                case 'Piston'
                    obj.order = 0;
                    obj.repetition = 0;
                    
                case 'Tip'
                    obj.order = 1;
                    obj.repetition = 1;
                    
                case 'Tilt'
                    obj.order = 1;
                    obj.repetition = -1;
                    
                case 'Defocus'
                    obj.order = 2;
                    obj.repetition = 0;
                    
                case 'Oblique astigmatism'
                    obj.order = 2;
                    obj.repetition = -2;
                    
                case 'Vertical astigmatism'
                    obj.order = 2;
                    obj.repetition = 2;
                    
                case 'Vertical coma'
                    obj.order = 3;
                    obj.repetition = -1;
                    
                case 'Horizontal coma'
                    obj.order = 3;
                    obj.repetition = 1;
                    
                case 'Vertical trefoil'
                    obj.order = 3;
                    obj.repetition = -3;
                    
                case 'Oblique trefoil'
                    obj.order = 3;
                    obj.repetition = 3;
                    
                case 'Primery spherical'
                    obj.order = 4;
                    obj.repetition = 0;
                         
                case 'Vertical secondary astigmatism'   
                    obj.order = 4;
                    obj.repetition = 2;
                    
                case 'Oblique secondary astigmatism'   
                    obj.order = 4;
                    obj.repetition = -2;
                    
                case 'Vertical quadrafoil'   
                    obj.order = 4;
                    obj.repetition = 4;
                    
                case 'Oblique quadrafoil'   
                    obj.order = 4;
                    obj.repetition = -4;
                    
                otherwise
                    error('Unknown aberration name');
            end
            obj.aberration_name = type;
        end
   end
   
   
   % private methods
   methods (Access = private)
       
        % -----------------------
        % set the aberration type util
        function obj = setAberrationTypeUtil(obj)
           switch obj.order
                case 0
                    obj.aberration_name = 'Piston';

                case 1
                    if obj.repetition == 1
                        obj.aberration_name = 'Tip';
                    elseif obj.repetition == -1
                        obj.aberration_name = 'Tilt';
                    end

                case 2
                    if obj.repetition == 0
                        obj.aberration_name = 'Defocus';
                    elseif obj.repetition == -2
                        obj.aberration_name = 'Oblique astigmatism';
                    elseif obj.repetition == 2
                        obj.aberration_name = 'Vertical astigmatism';
                    end

                case 3
                    if obj.repetition == -1
                        obj.aberration_name = 'Vertical coma';
                    elseif obj.repetition == 1
                        obj.aberration_name = 'Horizontal coma';
                    elseif obj.repetition == -3
                        obj.aberration_name = 'Vertical trefoil';
                    elseif obj.repetition == 3
                        obj.aberration_name = 'Oblique trefoil';
                    end

                case 4
                    if obj.repetition == 0
                        obj.aberration_name = 'Primery spherical';
                    elseif obj.repetition == 2
                        obj.aberration_name = 'Vertical secondary astigmatism';
                    elseif obj.repetition == -2
                        obj.aberration_name = 'Oblique secondary astigmatism';
                    elseif obj.repetition == 4
                        obj.aberration_name = 'Vertical quadrafoil';
                    elseif obj.repetition == -4
                        obj.aberration_name = 'Oblique quadrafoil';
                    end

                otherwise
                    obj.aberration_name = 'High order';
            end
       end
   end
end