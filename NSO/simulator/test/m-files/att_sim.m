clear all;
clc;


for k_att=1:30
        %Generating initial quaternion
        theta = rand*2*pi;
        phi = rand*2*pi;
        alpha = rand*2*pi;
        %Random e-vector
        [e(1) e(2) e(3)]=sph2cart(theta,phi,1);
        %Random quaternion used in the mdl mask of the NSO.
        q_init = [e(1)*sin(alpha/2) e(2)*sin(alpha/2) e(3)*sin(alpha/2) cos(alpha/2)];
       
        %Generating reference quaternion
        theta = rand*(pi/4);
        phi   = rand*(pi/4);
        alpha = rand*(pi/4);
        
        q_x =[sin(theta/2) 0 0 cos(theta/2)];
        q_y =[0 sin(phi/2) 0 cos(phi/2)];
        q_z = [0 0 sin(alpha/2) cos(alpha/2)];
                
        q_ref = qmulatt(q_init,q_x);
        q_ref = qmulatt(q_ref,q_y);
        q_ref = qmulatt(q_ref,q_z);
        
         OK=false;
         while(~OK)
         [OK q_ref]=checkatt(q_init,q_ref);
         end
        
        [T_att{k_att},X,Y_att{k_att}]=sim('test/attitude_control');
        
        q_att{k_att}=[theta phi alpha q_ref q_init];
k_att
end

save attitude_mcsim.mat T_att Y_att q_att;

