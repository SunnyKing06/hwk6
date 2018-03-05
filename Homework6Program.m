%PROGRAM FOR HWK6
%WRITTEN IN MATLAB; NOT in C
clear all
clc
close(gcf)
% Set up boundary Conditions and step increments
%also set up plate length
% M for x; N for y;
lx= pi;
ly=pi;
%set up nodes here
%MAKE SURE THAT IT IS EVEN
m=1000;
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
iteration=4000;
%defining x & y intervals

x=[dx:dx:lx-dx];
y=[dx:dy:ly-dy];

%Closed-form solution given in Assignment sheet 
% M an integer
M=2;
for p=1:length(x)
    for q=1:length(y)
      solU(p,q)=(ly-y(q)).*sin(M.*x(p)).*sinh(M.*y(q));
      solF(p,q)= -2.*M.*sin(M.*x(p)).*cosh(M.*y(q));
    end
end
surf(x,y,solU)
%given f in del2u=f

surf(x,y,solF);

%setting up u properly:
u(m,n)=0;
u_new(m,n)=0;
error(m,n)=0;
%to not compute too much in a loop, will do the following 
l= dx*dy*2*M;

for j=1:n
    for i=1:m
       %checking to see if at boundaries
        if i==1
           if j==1
             %will be at CORNER (0,0)  
                u(i,j)=(l*sin(M*i*dx)*cosh(M*j*dy) -u0-v0-u(i+1,j) -u(i,j+1))/-4;
           elseif j==n % will be at CORNER (0,ly)
                u(i,j)= (l*sin(M*i*dx)*cosh(M*j*dy) -u(i,j-1)-v0-u(i+1,j) -uL)/-4
           else %will be at LINE (0,y)
                u(i,j)= (l*sin(M*i*dx)*cosh(M*j*dy) -u(i,j-1)-v0-u(i+1,j) -u(i,j+1))/-4;
           end
        end
        
        if (j==1) && ((i~=1) && (i~=m)) % now will at LINE at (x !=0,0)
            u(i,j)= (l*sin(M*i*dx)*cosh(M*j*dy) -u0-u(i-1,j)-u(i+1,j) -u(i,j+1))/-4;
        end
         
        if i==m
            if j==n
                %will be at CORNER (lx,ly)
                 u(i,j)= (l*sin(M*i*dx)*cosh(M*j*dy) -u(i,j-1)-u(i-1,j)-uL -vL)/-4;
            elseif j==1
                % will be at CORNER (lx,0)
                u(i,j)= (l*sin(M*i*dx)*cosh(M*j*dy) -u0-u(i-1,j)-vL -u(i,j+1))/-4;
            else %now at LINE (lx,y)
                u(i,j)= (l*sin(M*i*dx)*cosh(M*j*dy) -u(i,j-1)-u(i-1,j)-vL -u(i,j+1))/-4;
            end
        end
        
         
         
         if (j==n) && ((i~=1) && (i~=m))  % will be at LINE (x,ly)
            u(i,j)= (l*sin(M*i*dx)*cosh(M*j*dy) -u(i,j-1)-u(i-1,j)-u(i+1,j) -uL)/-4;
         end
         
         %end of boundary check; else will be within space
        
        if ((i~=1) && (i~=m) && (j~=1) && (j~=n)) %checks to see if index does not depend on boundary    
            u(i,j)= (l*sin(M*i*dx)*cosh(M*j*dy) -u(i,j-1)-u(i-1,j)-u(i+1,j) -u(i,j+1))/-4;
        end
        %End of solving for u at boundaries & inner space area
        
     end
end
u_new=u;
%First guess done; will use gauss-seidel iteratively to get better
%answers
%again, to minimize multiplication within nested loops
l=-.25*dx*dy*-2*M;
        for k=1:iteration %repetition of each step
           for j=1:n
                for i=1:m
                    if (i==1) && (j==1) % lower  left corner check
                       u_new(i,j)= .25*(u0+u(i+1,j)+v0+u(i,j+1))+l*sin(M*i*dx)*cosh(M*j*dy);
                    
                    elseif (i==1) &&  (j==n) %upper left corner check 
                        u_new(i,j)= .25*(v0+u(i+1,j)+u_new(i,j-1)+uL)+l*sin(M*i*dx)*cosh(M*j*dy);
                    
                    elseif(i==m) &&(j==1) %lower right corner check
                        u_new(i,j)= .25*(u_new(i-1,j)+vL+u0+u(i,j+1))+l*sin(M*i*dx)*cosh(M*j*dy);
                        
                    elseif(i==m) &&(j==n) %upper right corner check
                        u_new(i,j)= .25*(u_new(i-1,j)+vL+u_new(i,j-1)+uL)+l*sin(M*i*dx)*cosh(M*j*dy);
                    
                    elseif (i==1) % x axis check
                        u_new(i,j)= .25*(v0+u(i+1,j)+u_new(i,j-1)+u(i,j+1))+l*sin(M*i*dx)*cosh(M*j*dy);
                    
                    elseif (i==m)  % x=lx check
                        u_new(i,j)= .25*(u_new(i-1,j)+u_new(i,j-1)+u(i,j+1))+l*sin(M*i*dx)*cosh(M*j*dy);
                    
                    elseif(j==1) % y axis check
                        u_new(i,j)= .25*(u_new(i-1,j)+u(i+1,j)+u(i,j+1))+l*sin(M*i*dx)*cosh(M*j*dy);
                    
                    elseif(j==m) %y=ly
                        u_new(i,j)= .25*(u_new(i-1,j)+u(i+1,j)+u_new(i,j-1))+l*sin(M*i*dx)*cosh(M*j*dy);
                    
                    else %all others wont depend on boundary
                        u_new(i,j)= .25*(u_new(i-1,j)+u(i+1,j)+u_new(i,j-1)+u(i,j+1))+l*sin(M*i*dx)*cosh(M*j*dy);
                    
                    end
                     u(i,j)=u_new(i,j);
                end
           end
           u=u_new;
        end
        
     u;
     solU;
     for i=1:m
         for j=1:n
             error(i,j)= abs((u(i,j)-solU(i,j))/solU(i,j))*100;
         end
     end
     error;
     naive_average_error= sum(sum(error))/(n*m)
     k=surf(x,y,solU)
     %makes true solution a little transparent so can see both surfaces
     % easier
     alpha(k,.4)
     hold on
     surf(x,y,u)
     
     %plotting random values of the approximated u and the Closed-form solU
%     s=m/2;
%     t=n/2;
%      rx= randi([1,m],1,s)
%     ry= randi([1,n],1,t)
%     
%     figure
%     plot(rx,u(rx,ry),'g') 
%     hold on
%     plot(rx, solU(rx,ry),'b')
