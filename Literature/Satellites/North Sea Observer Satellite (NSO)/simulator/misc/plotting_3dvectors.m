%B field plots
format long
global n_ZX n_XY n_YZ e Ltot Ltotend;% e_LL n_ZX_LL n_XY_LL n_YZ_LL;

size = 12;
e=(1/norm(e))*e;
Ltot=(1/norm(Ltot))*Ltot;
Ltotend=(1/norm(Ltotend))*Ltotend;
sun=[0.7852904564 0.0761112902 0.6144314206];
sun=(1/norm(sun))*sun;
nul = [0 0 0];


angle_ZX_LL=acos(dot(e,Ltot))
angle_eL=acos(dot(e,Ltotend))
angle_esun = acos(dot(e,sun))
angle_ex= acos(dot(e,[1 0 0]))
figure(1)
clf
arrow3(nul, [1 0 0], 'k'); nice3d
hold on
texts([1.1 0 0], 'x',size)
arrow3(nul, [0 1 0], 'k'); nice3d
texts([0 1.1 0], 'y',size)
arrow3(nul, [0 0 1], 'k'); nice3d
texts([0 0 1.1], 'z',size)

arrow3(nul, e, 'b'); nice3d
texts(e*1.1, 'e',size)
arrow3(nul, n_ZX, 'r'); nice3d
texts(n_ZX*1.1, 'n_ZX',size)
arrow3(nul, n_YZ, 'r'); nice3d
texts(n_YZ*1.1, 'n_YZ',size)
arrow3(nul, n_XY, 'r'); nice3d
texts(n_XY*1.1, 'n_XY',size)
arrow3(nul, Ltot, 'g'); nice3d
texts(Ltot*1.1, 'Ltot',size)
arrow3(nul, Ltotend, 'c'); nice3d
texts(Ltotend*1.1, 'Ltotend',size)

% figure(2)
% clf
% arrow3(nul, [1 0 0], 'k'); nice3d
% hold on
% texts([1.1 0 0], 'x',size)
% arrow3(nul, [0 1 0], 'k'); nice3d
% texts([0 1.1 0], 'y',size)
% arrow3(nul, [0 0 1], 'k'); nice3d
% texts([0 0 1.1], 'z',size)
% 
% arrow3(nul, e_LL, 'b'); nice3d
% texts(e_LL*1.1, 'Ltot x Ltotend',size)
% arrow3(nul, n_ZX_LL, 'r'); nice3d
% texts(n_ZX_LL*1.1, 'n_ZX',size)
% arrow3(nul, n_YZ_LL, 'r'); nice3d
% texts(n_YZ_LL*1.1, 'n_YZ',size)
% arrow3(nul, n_XY_LL, 'r'); nice3d
% texts(n_XY_LL*1.1, 'n_XY',size)
% arrow3(nul, Ltot, 'g'); nice3d
% texts(Ltot*1.1, 'Ltot',size)
% arrow3(nul, Ltotend, 'c'); nice3d
% texts(Ltotend*1.1, 'Ltotend',size)
