%% Define targets in the format {name, [lat,lon,alt]}
targets = {'sydnye', [-33.876433, 151.185346,0];
%            'wien', [48.206251, 16.373879, 542];
%            'sao_paulo' [-2.000000, -44.600000,10];
%            'bremen', [53.059346, 8.807349,11];
%            'aalborg HQ', [57.013835, 9.987509,15];
%            'Madrid', [40.415246, -3.707029, 667];
%            'Lisboa', [38.718298, -9.142485, 200];
%            'Paris', [48.840797, 2.351053, 35];
%            'Tuhle Basen' [76.5322, -68.0731, 500];
           'Paderborn' [51.5334, 8.8725, 10]
           };
       
save('target_list.mat','targets');