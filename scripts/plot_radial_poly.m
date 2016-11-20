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
%   Calculates and plots radial polynomials of different orders 
%
% Revision history:
%   7/13/2016 - Created from existing code.
%   11/19/2016 - Minor cosmatic enhancements.
%
% -------------------------------------------------------------------------
clc; close all;
% user input segment

n = 0:4;            % range of the desired orders for radial polynomials

rho = 0:0.01:1;     % range of radial distance 

% end user input segment
% -------------------------------------------------------------------------
frm = {'-r', '--k', '-g', '--g', '-b', '--b', '-m', '--m', '-.m'};
cnt = 1;
% initialize m
data = [];
% dummy plot (place holder for the first entry of the legend)
plot(0,0, 'w'); hold on; grid on;

for i=1:length(n)
    for m=0:n(i)
        if ~mod(n(i)-abs(m), 2)
            myZernike = zernike(n(i),m);
            data = [data; n(i),m];         %#ok
            R = myZernike.getRadialPoly(rho);
            plot(rho, R, frm{cnt}, 'Linewidth', 1);
            cnt = cnt + 1;
        end
    end
end
ylabel('Radial Polynomial', 'fontsize', 12);
xlabel('\rho', 'fontsize', 12);
legend(['n, m'; cellstr(num2str(data))]);
