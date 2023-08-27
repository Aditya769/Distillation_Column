clc
clear
close
options=optimset('display','off');



x1  = (0   :0.01:1);
x2 = 1-x1 ;
%disp(x2);
y1 = (length(x1));
A21 = -0.5899;
A12  = -0.8643;
T = 337.7; % in Kelvin
den = 1268; % density -----> in  Kg/m^3
P  = 100000;% in Pascal; 
% Calculating the fugacity coefficient 

% Acetone 
Tc =  781 ; % k      
Pc =   469000; % in Pascal Unit ;
R = 8.314;
a = 27*(R^2)(Tc^2)/(64(Pc^2));
b =  R*Tc/(8*Pc);

%Chloroform 
Tc2 =  536.4 ; % k      
Pc2 =  5400000;% in Pascal;
a2 = 27*(R^2)(Tc2^2)/(64(Pc2^2));
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
% Antonie 's Constants.........
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

% We use Van- Laar Model to compute  activity coefficient;
for i =  1:length(x1)
gama1 = exp(A12*(A21*(x2(i))/(A12*x1(i)+A21*x2(i)))^2);
gama2 = exp(A21*(A12*(x1(i))/(A12*x1(i)+A21*x2(i)))^2);
y1(i)  = 1/((RATIO*gama2*x2(i))/(k*gama1*x1(i))+1);
%disp(y1(i));
end
t =  (RATIO*gama2)/(gama1*k);
% Taking the input from user.............

% You will need the value of reflux ratio, boil-up ratio, 
% bottom and top product composition to solve the above problems.
% xd, xf, xw ,q ;
disp("Please enter all the composition (mole-fraction) of more volatile components in the following questions");
x_top= input("Enter the composition of top section");
x_feed= input("Enter the composition of feed ");
x_bottom =input("Enter the composition of bottom");
q = input("Enter the value of q ");




% feed Line
Q= (q/(q-1));
C2 = x_feed/(q-1);

if q>1
    xq2=1;
    yq2=Q*xq2-C2;
elseif q==1
    xq2=x_feed;
    yq2=1;
elseif q<1 && q>0
    xq2=0;
    yq2=Q*xq2-C2;
elseif q==0
    xq2=0;
    yq2=x_feed;
else
    xq2=0;
    yq2=Q*xq2-C2;
end
    
y_feed  = @(x) Q*x - C2;

% Calculating the minimum Reflux Ratio
%Equlibrium equation =>  y = xk/(xk+1-x)
eqlbrm_eq=  @(x) x/(t*(1-x)+x);

if q==1
    x_common=x_feed;
    y_common=eqlbrm_eq(x_common);
elseif q==0
    y_common=x_feed;
    x_common= fsolve(@(x) y_common- x/(t*(1-x)+x),x_feed,options);
else
    x_common=fsolve(@(x)  x/(t*(1-x)+x) -(Q*x-C2),0.5,options);
    y_common=eqlbrm_eq(x_common);
end
R_min_slope=(x_top-y_common)/(x_top-x_common);

R_min_intercept = x_top - R_min_slope*x_top;
R_min = x_top/R_min_intercept -1;
reflux_ratio= 1.2*R_min;
disp("Minimum Reflux ratio comes out to be:")
disp(R_min);

% Enriching Section
R = (reflux_ratio/(reflux_ratio + 1));
C1 =x_top/(reflux_ratio+1);



% %Finding the intersection point between Enriching line and Equilibria
x_intersect =  fsolve(@(x) x/(t*(1-x)+x)-(R*x+C1),x_top,options);
y_intersect = eqlbrm_eq(x_intersect);


  


R = (reflux_ratio/(reflux_ratio + 1));
C1 =x_top/(reflux_ratio+1);

y_top =@(x) R*x + C1;



%Bottom line 

% We use the 2 points  to find out equation : 
% (x_bottom,x_bottom) 
% and other point of 
%intersection of enriching line and feed line (a,b)  ..
if q==1
    a= x_feed;
    b=  R*a + C1;
else 

    a =  (-C2-C1)/(R-Q);
    b = Q*a - C2;
end    
slope = (b-x_bottom)/(a-x_bottom);
y_bottom = @(x) slope*(x-x_bottom) + x_bottom;

if a>x_intersect
    x_intersect=a;
end

% %Finding the intersection point between Stripping line and Equilibria
x_intersect_bot =  fsolve(@(x) x/(t*(1-x)+x)-(slope*(x-x_bottom) + x_bottom),0.4,options);
y_intersect_bot = eqlbrm_eq(x_intersect_bot);

if b<y_intersect_bot
    y_intersect_bot =b;

end

% PLOTTING OF Equilibrium Equation :
plot(x1,y1,'-k');
hold on

% 45 degree line 
plot(x1,x1,'-k');
xlabel('x'),ylabel('y')
set(gca,'Xlim',[0 1])
set(gca,'Ylim',[0 1])
hold on
grid on 

% Drawing top section line 

plot([x_intersect x_top],[y_intersect x_top],'-k','linewidth',0.5);
hold on 

%Plotting feed line

plot([x_feed xq2],[x_feed yq2],'-k','linewidth',0.5);
hold on 
% 
%Plotting bottom line 
plot([x_bottom a],[x_bottom b],'-k');

% legend('Equilibrium Curve','straight line','enriching line','feed line','stripping line','location','northwest');

% Top section stages
% 
x_eqn  = @(y) t*y/(1-y+t*y) ;
topsection = @(x) R*x + C1;
legend('off')

x_top_1 = x_top;
y_top_1 = x_top;
i=0;
while x_top_1>x_intersect
    y_top_2 = y_top_1;
    x_top_2 = x_eqn(y_top_2);
    x_top_3 = x_top_2;
    y_top_3 = topsection(x_top_3);
    
    plot([x_top_1 x_top_2],[y_top_1 y_top_2],'r' );
    
    plot([x_top_2 x_top_3],[y_top_2 y_top_3],'r' );
    
    x_plot_top = x_top_1;
    x_top_1=x_top_3;
    y_top_1=y_top_3;
    i=i+1;

end
Stages_enriching_section= i-(x_top_2-x_intersect)/(x_top_2-x_plot_top); 
disp("Stages in enriching Section");
disp(Stages_enriching_section);
% 
% STRIPPING SECTION

% % % % % x_bot = @(y) (y-x_bottom)/slope + x_bottom;
% % % % % 
% % % % % x_bottom_1 = x_bottom;
% % % % % y_bottom_1 = x_bottom;
% % % % % j = 0;
% % % % % 
% % % % % while y_bottom_1<b
% % % % %     x_bottom_2 = x_bottom_1;
% % % % %     y_bottom_2 = eqlbrm_eq(x_bottom_2);
% % % % %     y_bottom_3 = y_bottom_2;
% % % % %     x_bottom_3 = x_bot(y_bottom_3);
% % % % %     plot([x_bottom_1 x_bottom_2],[y_bottom_1 y_bottom_2], 'b');
% % % % %     plot([x_bottom_2 x_bottom_3],[y_bottom_2 y_bottom_3], 'b');
% % % % %     y_bottom_plot = y_bottom_1;
% % % % %     x_bottom_1 = x_bottom_3;
% % % % %     y_bottom_1 = y_bottom_3;
% % % % %     j = j+1;   
% % % % % end
% % % % % 
% % % % % Stages_stripping_section = j-(y_bottom_2-b)/(y_bottom_2-y_bottom_plot);
% % % % % disp("Total Number of Stages for  the given compositions are:");
% % % % % disp(Stages_enriching_section+Stages_stripping_section);