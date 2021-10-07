function [fx] = Obj(x)
% Phase 2 : Defining objective function
x1=x(:,1);
x2=x(:,2);
fx=x1.^2-x1.*x2+x2.^2+2.*x1+4.*x2+3;
end

