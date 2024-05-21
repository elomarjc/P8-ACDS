function [bw_robust bw_simple] = calbw(I_in)
% Calculaes BW of the closed loop systems with the defined inertia matrix (I_s).
% Runs all the design files needed to do so.

%Robust control
ang_robust_control;
att_control_with_integrator;
bw_robust = getbw(D_att_c,D_ang_c,I_in);

%Simpel control
basic_ang_control;
basic_attitude_control;
bw_simple = getbw(D_att_c,D_ang_c,I_in);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bw = getbw(K_att,K_ang,I_ss)
% Calculaes BW of a closed loop system in the from of the structure in the
% report
allowed_bw = 1.256;

operating_point_constants;

A_sat=[inv(I_ss)*Sh];
B_sat=-inv(I_ss);
C_sat=eye(3);
D_sat=zeros(3);
SYS_ang=ss(A_sat,B_sat,C_sat,D_sat);
CL_ang=feedback(K_ang*SYS_ang,eye(3),+1);


A_att = [zeros(3) 0.5*eye(3); zeros(3) CL_ang.a];
B_att = [zeros(3); CL_ang.b];
C_att = [eye(3) zeros(3);zeros(3) CL_ang.c];
D_att = zeros(6,3);

if(size(K_att,2)==9)
    %integral action
    A_att = [A_att zeros(6,3); eye(3) zeros(3,6)];
    B_att = [B_att; zeros(3)];
    C_att = eye(9);
    D_att = zeros(9,3);
end

SYS_att=ss(A_att,B_att,C_att,D_att);
CL_att=feedback(K_att*SYS_att,eye(3),+1);

%[Wn,Z,P] = damp(CL_att);
%resp = FREQRESP(CL_att,[1.0:0.01:2]);

resp = FREQRESP(CL_att,allowed_bw);
db =[];
for it=1:length(resp)
    db = [db 20*log10(norm(resp(it,it)))];
end
db=max(db);
step=0.001;
idx=1;
if db<=-3
    while db<=-3
        resp = FREQRESP(CL_att,allowed_bw-step*idx);
        db = [];
        for it=1:length(resp)
            db = [db 20*log10(norm(resp(it,it)))];
        end
        db=max(db);
        idx = idx+1;
    end
    bw = allowed_bw-step*(idx-1);
else
    while db>-3
        resp = FREQRESP(CL_att,allowed_bw+step*idx);
        db = [];
        for it=1:length(resp)
            db = [db 20*log10(norm(resp(it,it)))];
        end
        db=max(db);
        idx = idx+1;
    end
    bw = allowed_bw+step*(idx-1);
end
bw = [bw bw/(2*pi)];

%bw = [Wn,Z,P];
%bw = max(Wn)/(2*pi);



