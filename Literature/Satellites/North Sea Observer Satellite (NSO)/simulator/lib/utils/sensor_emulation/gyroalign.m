function return_value = gyroalign(angle)


%gyro1
        err=angle(1)*(pi/180);
        n=[1;0;0];	
		C2=[cos(err) 0 -sin(err);0 1 0; sin(err) 0 cos(err)];
		C3=[cos(err) sin(err) 0; -sin(err) cos(err) 0; 0 0 1];
		Cres=C2*C3;
		gyro1=Cres*n;

%gyro2
        err=angle(2)*(pi/180);
        n=[0;1;0];
		C1=[1 0 0;0 cos(err) sin(err); 0 -sin(err) cos(err)];
		C3=[cos(err) sin(err) 0; -sin(err) cos(err) 0; 0 0 1];
		Cres=C1*C3;
		gyro2=Cres*n;

%gyro3
        err=angle(3)*(pi/180);
        n=[0;0;1];
		C1=[1 0 0;0 cos(err) sin(err); 0 -sin(err) cos(err)];
		C2=[cos(err) 0 -sin(err);0 1 0; sin(err) 0 cos(err)];
		Cres=C1*C2;
		gyro3=Cres*n;

%gyro4    
        err=angle(4)*(pi/180);
        n=[1;0;0];	
		C2=[cos(err) 0 -sin(err);0 1 0; sin(err) 0 cos(err)];
		C3=[cos(err) sin(err) 0; -sin(err) cos(err) 0; 0 0 1];
		Cres=C2*C3;
		gyro4=Cres*n;

%gyro5
        err=angle(5)*(pi/180);
        n=[0;1;0];
		C1=[1 0 0;0 cos(err) sin(err); 0 -sin(err) cos(err)];
		C3=[cos(err) sin(err) 0; -sin(err) cos(err) 0; 0 0 1];
		Cres=C1*C3;
		gyro5=Cres*n;

%gyro6
        err=angle(6)*(pi/180);
        n=[0;0;1];
		C1=[1 0 0;0 cos(err) sin(err); 0 -sin(err) cos(err)];
		C2=[cos(err) 0 -sin(err);0 1 0; sin(err) 0 cos(err)];
		Cres=C1*C2;
		gyro6=Cres*n;

        return_value=[gyro1(1) gyro2(2) gyro3(3) gyro4(1) gyro5(2) gyro6(3)];