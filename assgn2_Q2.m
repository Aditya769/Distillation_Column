clc
clear
close
options=optimset('display','off');

x = (0.01:0.01:1);
y = (length(x)); 
P_acetone_sat  = 0.30;%0.30;%230.6;%345;  % (For temperature range --(20 -25 degree Celcius)
P_chloroform_sat =  0.25;%0.25;%196.6; %200; 
k = P_acetone_sat/P_chloroform_sat;
for i = 1:length(x)
eqn = @(y) y-x(i)*y-k*x(i)+k*x(i)*y;
y(i) = fsolve(eqn,0.5);
disp(y(i));
end


% Taking the input from user.............


disp("Please enter all the composition (mole-fraction) of more volatile components in the following questions");
x_top= input("Enter the composition of top section");
x_feed= input("Enter the composition of feed ");
x_bottom =input("Enter the composition of bottom");
q = input("Enter the value of q ");
% reflux_ratio = 2.3;
% x_top=0.9;
% x_feed=0.7;
% x_bottom =0.05;
% q =1;





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
% Equlibrium equation =>  y = xk/(xk+1-x)
eqlbrm_eq=  @(x) x*k/(x*k +1-x);
if q==1
    x_common=x_feed;
    y_common=eqlbrm_eq(x_common);
elseif q==0
    y_common=x_feed;
    x_common= fsolve(@(x) y_common-(x*k/(x*k +1-x)),x_feed,options);
else
    x_common=fsolve(@(x)(x*k/(x*k +1-x)) -(Q*x-C2),0.5,options);
    y_common=eqlbrm_eq(x_common);
end
R_min_slope=(x_top-y_common)/(x_top-x_common);
R_min_intercept = x_top - R_min_slope*x_top;
R_min = x_top/R_min_intercept -1;
reflux_ratio= 1.3*R_min;


% Enriching Section
R = (reflux_ratio/(reflux_ratio + 1));
C1 =x_top/(reflux_ratio+1);

y_top = R*x + C1;
% Bottom line 

% We use the 2 points  to find out equation : 
% (x_bottom,x_bottom) 
% and other point of 
% intersection of enriching line and feed line (a,b)  ..
if q==1
    a= x_feed;
    b=  R*a + C1;
else 

    a =  (-C2-C1)/(R-Q);
    b = Q*a - C2;
end    
slope = (b-x_bottom)/(a-x_bottom);
y_bottom =  slope*(x-x_bottom) + x_bottom;

% 
nexttile;
% PLOTTING OF Equilibrium Equation :
plot(x,y,'-r');
hold on

% 45 degree line 
plot(x,x,'-g');
xlabel('x'),ylabel('y')
set(gca,'Xlim',[0 1])
set(gca,'Ylim',[0 1])
hold on
grid on 

% Drawing top section line 

plot([a x_top],[b x_top],'-b','linewidth',0.5);
hold on 

% Plotting feed line

plot([x_feed xq2],[x_feed yq2],'-k','linewidth',0.5);
hold on 

% Plotting bottom line 
plot([x_bottom a],[x_bottom b],'-m');

% legend('Equilibrium Curve','straight line','enriching line','feed line','stripping line','location','northwest');
legend('off');
% Top section stages
x_eqn  = @(y) y/(y+k-k*y);
topsection = @(x) R*x + C1;
% legend('off')
% 


x_top_1 = x_top;
y_top_1 = x_top;
i=0;

while x_top_1>a
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
    
 Stages_enriching_section= i-(x_top_2-a)/(x_top_2-x_plot_top); 

% %STRIPPING SECTION

y_eqn_ = @(x) x*k/(k*x+1-x);
x_bot = @(y) (y-x_bottom)/slope + x_bottom;

x_bottom_1 = x_bottom;
y_bottom_1 = x_bottom;
j = 0;
while y_bottom_1<b
    x_bottom_2 = x_bottom_1;
    y_bottom_2 = y_eqn_(x_bottom_2);
    y_bottom_3 = y_bottom_2;
    x_bottom_3 = x_bot(y_bottom_3);
    plot([x_bottom_1 x_bottom_2],[y_bottom_1 y_bottom_2], 'b');
    plot([x_bottom_2 x_bottom_3],[y_bottom_2 y_bottom_3], 'b');
    y_bottom_plot = y_bottom_1;
    x_bottom_1 = x_bottom_3;
    y_bottom_1 = y_bottom_3;
    j = j+1;   
end
Stages_stripping_section = j-(y_bottom_2-b)/(y_bottom_2-y_bottom_plot);
disp("Total Number of Stages for  the given compositions are:");
disp(Stages_enriching_section+Stages_stripping_section);
hold off;

nexttile;
plot(x,y,'-r');
hold on;
plot(x,x,'-g');
xlabel('x'),ylabel('y')
set(gca,'Xlim',[0 1])
set(gca,'Ylim',[0 1])
hold on
grid on
legend('off');
% For Calculating minimum number of trays:
%reflux ratio = INFINITE;
topsection_min_trays = @(x) x;



x_top_1 = x_top;
y_top_1 = x_top;
p=0;

while x_top_1>x_bottom
    y_top_2 = y_top_1;
    x_top_2 = x_eqn(y_top_2);
    x_top_3 = x_top_2;
    y_top_3 = topsection_min_trays(x_top_3);
    
    plot([x_top_1 x_top_2],[y_top_1 y_top_2],'r' );
    
    plot([x_top_2 x_top_3],[y_top_2 y_top_3],'r' );
    
    x_plot_top_1 = x_top_1;
    x_top_1=x_top_3;
    y_top_1=y_top_3;
    p=p+1;

end
 disp("MINIMUM NUMBER OF TRAYS ARE :") ;
 disp(p-(x_top_2-x_bottom)/(x_top_2-x_plot_top_1)); 
 disp("Minimum Reflux Ratio comes out to be :")
 disp(R_min);