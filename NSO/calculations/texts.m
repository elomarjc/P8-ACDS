function texts(P,s,fontsize) 
% texts(P,'text') 
% 
% plot 'text' at point P in graph, P is vector of length 2 or 3 
P = double(P); 
if length(P)==2 
  H=text(P(1),P(2),s); 
elseif length(P)==3 
  H=text(P(1),P(2),P(3),s); 
else 
  error('P must be vector of length 2 or 3') 
end 
 
set(H, 'FontSize', fontsize);