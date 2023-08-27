x1  = (0 :0.01:1);
x2 = 1-x1;
%disp(x2);
y1 = (length(x1));
R = 8.314;
A21 = -0.5899;
A12  = -0.8643;
v1 = 74.166; % molar volumes of acetone in cm3/mol
v2 = 80.235; % molar volumes of chloroform in cm3/mol
T = 337.7; % in Kelvin
den = 1268; % density -----> in  Kg/m^3
P  = 100000;% in Pascal; 
% calculation of a12 and a21 for wilson model where
% a12 = (v2/v1)*exp(-A12/RT), a21 = (v1/v2)*exp(-A21/RT)
a12 = (v2/v1)*exp(-A12/R*T);
a21 = (v1/v2)*exp(-A21/R*T);
% Calculating the fugacity coefficient 

% Acetone 
Tc =  781 ; % k      
Pc =   469000; % in Pascal Unit ;
a = 27*(R^2)*(Tc^2)/(64*(Pc^2));
b =  R*Tc/(8*Pc);

%Chloroform 
Tc2 =  536.4 ; % k      
Pc2 =  5400000;% in Pascal;
a2 = 27*(R^2)*(Tc2^2)/(64*(Pc2^2));
b2 =  R*Tc2/(8*Pc2);


q1 = a/(b*R*T);
Z1 =  1/(1-den*b) - (a*den)/(R*T);
B1 = b*P/(R*T);

fug_Coeff_1 =  exp( Z1-1 -log(abs(Z1 - B1)) - q1*B1/Z1 );
disp(fug_Coeff_1);

q2= a2/(b2*R*T);
Z2 =  1/(1-den*b2) - (a2*den)/(R*T);
B2 = b2*P/(R*T);

fug_Coeff_2 =  exp( Z2-1 -log(Z2 - B2) - q2*B2/Z2 );
%disp(fug_Coeff_2);

RATIO = fug_Coeff_1/ fug_Coeff_2;
% disp(RATIO);
AA= 4.42448;
AB= 1312.253;
AC =32.445;

P1_sat =  (10^(AA- (AB/(T+AC))));  % bar 

AA2= 4.20772;
AB2= 1233.129;
AC2 =-40.953;

P2_sat =  (10^(AA2- (AB2/(T+AC2)))) ;

k = P1_sat/P2_sat;
disp(k);

%y1  = 1/((RATIO*gama2*x2)/(k*gama1*x1)+1);

% We use wilson Model to compute  activity coefficient;
for i =  1:length(x1)
gama1 = exp(-log((1-x2(i))+a12*x2(i))+x2(i)*((a12/((1-x2(i))+a12*x2(i)))-(a21/(x2(i)+a21*(1-x2(i))))));
gama2 = exp(-log(x2(i)+a21*(1-x2(i)))+(1-x2(i))*((a21/((x2(i))+a12*(1-x2(i))))-(a12/((1-x2(i))+a12*(x2(i))))));
y1(i)  = 1/((RATIO*gama2*x2(i))/(k*gama1*(1-x2(i)))+1);
disp(y1(i));
end

v = x1;
plot(x1,y1);
hold on 
plot(x1,v);