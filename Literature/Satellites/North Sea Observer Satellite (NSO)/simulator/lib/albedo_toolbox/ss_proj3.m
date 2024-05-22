% SS_PROJ Project all cell albedo contributions onto the solar cell normal
% in ECEF. This value is equivalent to the total perpendicular irradiance
% reaching the solar cell.
%
% P = ss_proj(re,nx,ny,a)
%
% where a is an albedo matrix, nx and ny are the cell normal to the x- and y-axis in ECEF frame, and
% re is the nadir in ECEF frame.
%
% $Id: ss_proj3.m,v 1.2 2006/04/18 14:27:19 mnkr02 Exp $

function P = ss_proj(re,nx,ny,a);
n = [nx,ny,cross(nx,ny)];
CONST.EMR = 6371.01e3;
CONST.AM0 = 1366.9;

P = [0 0 0 0 0 0];

% function [theta,phi] = idx2rad(i,j,sy,sx);
%
% CONST.d2r = pi/180;
%
% dx = 2*pi/sx;
% dy = pi/sy;
%
% phi = pi-dy/2-(i-1)*dy;
% theta = (j-1)*dx-pi+dx/2;
%
persistent grid;
[sy sx] = size(a);
if isempty(grid)
    grid = zeros(3,sy,sx);
    for r=1:sy
        for t=1:sx
            [grid_theta grid_phi]= idx2rad(r,t,sy,sx);
            [grid(1,r,t) grid(2,r,t) grid(3,r,t)] = sph2cart(grid_theta,pi/2-grid_phi,CONST.EMR);
        end
    end
end
for j=1:sx
    if norm(a(:,j))>0
        for i=1:sy
            if a(i,j) > 0
                % Grid vector in ECEF
                %[grid_theta grid_phi] = idx2rad(i,j,sy,sx);
                %[grid(1) grid(2) grid(3)] = sph2cart(grid_theta,pi/2-grid_phi,CONST.EMR);
                % Cosine of the angle from solar cell normal to grid LOS vector
                % (re+grid)
                grid_vec = (re+grid(:,i,j))/norm(re+grid(:,i,j));
                cosphi = n'*grid_vec;
                for k=1:3
                    if cosphi(k) > 0
                        P(k) = P(k) + a(i,j)*cosphi(k);
                    else
                        P(k+3) = P(k+3) + a(i,j)*-cosphi(k);
                    end
                end
            end
        end
    end
end