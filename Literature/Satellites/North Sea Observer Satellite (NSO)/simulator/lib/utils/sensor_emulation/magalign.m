function res = magalign(angle)
err=angle*(pi/180);
		C1=[1 0 0;0 cos(err) sin(err); 0 -sin(err) cos(err)];
		C2=[cos(err) 0 -sin(err);0 1 0; sin(err) 0 cos(err)];
		C3=[cos(err) sin(err) 0; -sin(err) cos(err) 0; 0 0 1];
		res=C1*C2*C3;