%y =(10);
%disp(y);
%x = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
x = (0.01:0.01:1);
y = (length(x)); 
P_acetone_sat  = 30;  % (For temperature range --(20 -25 degree Celcius)
P_chloroform_sat =  26; 
k = P_acetone_sat/P_chloroform_sat;
for i = 1:length(x)
eqn = @(y) y-x(i)*y-k*x(i)+k*x(i)*y;
y(i) = fsolve(eqn,0.5);
disp(y(i));
end
plot(x,y);
hold on
plot(x,x);