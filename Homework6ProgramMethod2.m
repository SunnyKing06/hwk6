%PROGRAM FOR HWK6
%WRITTEN IN MATLAB; NOT in C
clear all
clc
close(gcf)
% Set up boundary Conditions and step increments
%also set up plate length
% M for x; N for y;
lx= pi;
ly= pi;

%set up nodes here
%MAKE SURE THAT IT IS EVEN
m=100;
n=m;
dx= lx/(m+1);
dy= ly/(n+1);
%Boundary Conditions Defined here
u0=0;
uL=0;
v0=0;
vL=0;
%determining how many iterations to do to solver using 
%Gauss Seidel
iteration=600;
%defining x & y intervals

x=[0:dx:lx];
y=[0:dy:ly];

%Closed-form solution given in Assignment sheet 
% M an integer
M=1;
for p=1:length(x)
    for q=1:length(y)
      solU(p,q)=(ly-y(q)).*sin(M.*x(p)).*sinh(M.*y(q));
      solF(p,q)= -2.*M.*sin(M.*x(p)).*cosh(M.*y(q));
    end
end
surf(x,y,solU);
%given f in del2u=f

surf(x,y,solF);

%setting up u properly:
u(m+2,n+2)=0; % will bc's on the 'plate'
u(:,1)=u0;
u(:,m+2)=uL;
u(1,:)=v0;
u(n+2,:)=vL;

L1_error=0;
error=0;
u_new(m+2,n+2)=0;

%This nested loop solves for the first u values
%to not compute too much in a loop, will do the following 
l1= -dx*dy*2*M;
for s=1:iteration
    for j=2:n+1
        for i=2:m+1
       
       u(i,j)= (l1*sin(M*(i-1)*dx)*cosh(M*(j-1)*dy) -u_new(i,j-1)-u_new(i-1,j)-u(i+1,j) -u(i,j+1))/-4;
        u_new(i,j)=u(i,j); % new value assigned as in gauss-seidel
        end %end of i loop
        
    end % end of j loop
   u=u_new;      %iteration happens here
end

%Calculating the L1 error
     for i=2:m+1
         for j=2:n+2
             error= abs(u(i,j)-solU(i,j));
           L1_error=L1_error+error;
         end
     end
     L1_error= L1_error/(m*n)
     k=surf(x,y,solU)
     %makes true solution a little transparent so can see both surfaces
     % easier
     alpha(k,.4)
     hold on
     surf(x,y,u)
     
     %plotting random values of the approximated u and the Closed-form solU
%     
%      rx= randi([1,m],1,4)
%     ry= randi([1,n],1,4)
%     
%     figure
%     plot(rx,u(rx,ry),'g') 
%     hold on
%     plot(rx, solU(rx,ry),'b')
