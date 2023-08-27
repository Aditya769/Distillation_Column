
y =(10);
x = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
k = 7.99/1.64;
for i = 1:10
eqn = @(y) y-x(i)*y-k*x(i)+k*x(i)*y;
y(i) = fsolve(eqn,0.5);
%disp(y(i));
end
plot(x,y);