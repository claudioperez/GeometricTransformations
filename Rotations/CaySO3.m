function R = CaySO3(x)

theta=norm(x);
c1=4.0/(4.0+theta^2); 
c2=0.5*c1; % Cayley map coefficients

R = [1.0+c2*(x(1)^2-theta^2) c2*x(1)*x(2)-c1*x(3)  c2*x(1)*x(3)+c1*x(2);
    c2*x(1)*x(2)+c1*x(3)   1.0+c2*(x(2)^2-theta^2) c2*x(2)*x(3)-c1*x(1);
    c2*x(1)*x(3)-c1*x(2)     c2*x(2)*x(3)+c1*x(1) 1.0+c2*(x(3)^2-theta^2)]; 

