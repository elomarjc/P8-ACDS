function nice3d
% make 3d plots nicer by using perspective, showing bounding box, making
% axis units equal, showing grid, allow rotation by mouse dragging
xlabel('x'); ylabel('y'); zlabel('z')
camproj('perspective')    % undo by   camproj('orthographic')
box on                    % undo by   box off
axis tight                %
axis equal                % undo by   axis normal 
grid on                   % undo by   grid off

axis vis3d

%rotate3d on
cameratoolbar('SetMode','orbit');


