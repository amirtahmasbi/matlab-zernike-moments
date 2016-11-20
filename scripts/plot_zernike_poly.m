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
%   Calculates and plots Zernike polynomials of different orders 
%
% Revision history:
%   11/20/2016 - Created from existing code.
%
% -------------------------------------------------------------------------
clc; close all; clear all; %#ok
% user input segment

n   = 0:4;            % range of the desired orders for radial polynomials
roi = [31, 31];
% -----------------

for i=1:length(n)
    for m=-n(i):n(i)
        if ~mod(n(i)-abs(m), 2)
            myZernike = zernike(n(i),m);
            [moment, pupil] = myZernike.getZernikeMoment(zeros(roi));
            figure(1);subplot(length(n),2*length(n), n(i)*2*length(n)+m+length(n));
            mesh(pupil); ...colormap default;
            type = myZernike.getAberrationType();
            title(type);
            axis off; axis equal;
        end
    end
end
