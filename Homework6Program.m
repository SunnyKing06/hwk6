%PROGRAM FOR HWK6
%WRITTEN IN MATLAB
clear all
clc
% Set up boundary Conditions and step increments
%also set up plate length
% M for x; N for y;
lx= pi;
ly=pi;
m=10;
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
iteration=10;
%defining x & y intervals

x=[0:dx:lx];
y=[0:dy:ly];

%Closed-form solution given in Assignment sheet 
% M an integer
M=1;
for p=1:length(x)
    for q=1:length(y)
      solu(p,q)=(ly-y(q)).*sin(M.*x(p)).*sinh(M.*y(q));
      solf(p,q)= -2.*M.*sin(M.*x(p)).*cosh(M.*y(q));
    end
end
surf(x,y,solu);
%given f in del2u=f

%surf(x,y,solf);

%setting up u properly:
u(m,n)=0;

u_new(m,n)=0;

for j=1:n
    for i=1:m
       %checking to see if at boundaries
        if i==1
           if j==1
             %will be at CORNER (0,0)  
                u(i,j)=(dx*dy*2*M*sin(M*i*dx)*cosh(M*j*dy) -u0-v0-u(i+1,j) -u(i,j+1))/-4;
           elseif j==n % will be at CORNER (0,ly)
                u(i,j)= (dx*dy*2*M*sin(M*i*dx)*cosh(M*j*dy) -u(i,j-1)-v0-u(i+1,j) -uL)/-4
           else %will be at LINE (0,y)
                u(i,j)= (dx*dy*2*M*sin(M*i*dx)*cosh(M*j*dy) -u(i,j-1)-v0-u(i+1,j) -u(i,j+1))/-4;
           end
        end
        
        if (j==1) && ((i~=1) && (i~=m)) % now will at LINE at (x !=0,0)
            u(i,j)= (dx*dy*2*M*sin(M*i*dx)*cosh(M*j*dy) -u0-u(i-1,j)-u(i+1,j) -u(i,j+1))/-4;
        end
         
        if i==m
            if j==n
                %will be at CORNER (lx,ly)
                 u(i,j)= (dx*dy*2*M*sin(M*i*dx)*cosh(M*j*dy) -u(i,j-1)-u(i-1,j)-uL -vL)/-4;
            elseif j==1
                % will be at CORNER (lx,0)
                u(i,j)= (dx*dy*2*M*sin(M*i*dx)*cosh(M*j*dy) -u0-u(i-1,j)-vL -u(i,j+1))/-4;
            else %now at LINE (lx,y)
                u(i,j)= (dx*dy*2*M*sin(M*i*dx)*cosh(M*j*dy) -u(i,j-1)-u(i-1,j)-vL -u(i,j+1))/-4;
            end
        end
        
         
         
         if (j==n) && ((i~=1) && (i~=m))  % will be at LINE (x,ly)
            u(i,j)= (dx*dy*2*M*sin(M*i*dx)*cosh(M*j*dy) -u(i,j-1)-u(i-1,j)-u(i+1,j) -uL)/-4;
         end
         
         %end of boundary check; else will be within space
        
        if ((i~=1) && (i~=m) && (j~=1) && (j~=n)) %checks to see if index does not depend on boundary    
            u(i,j)= (dx*dy*2*M*sin(M*i*dx)*cosh(M*j*dy) -u(i,j-1)-u(i-1,j)-u(i+1,j) -u(i,j+1))/-4;
        end
        %End of solving for u at boundaries & inner space area
        
     end
end
u_new=u;
%First guess done; will use gauss-seidel iteratively to get better
%answers 
        for k=1:iteration %repetition of each step
           for j=1:n
                for i=1:m
                    if (i==1) && (j==1) % lower  left corner check
                       u_new(i,j)= .25*(u0+u(i+1,j)+v0+u(i,j+1))+-1*.25*dx*dy*-2*M*sin(M*i*dx)*cosh(M*j*dy);
                    
                    elseif (i==1) &&(j==n) %upper left corner check 
                        u_new(i,j)= .25*(v0+u(i+1,j)+u_new(i,j-1)+uL)+-1*.25*dx*dy*-2*M*sin(M*i*dx)*cosh(M*j*dy);
                    
                    elseif(i==m) &&(j==1) %lower right corner check
                        u_new(i,j)= .25*(u_new(i-1,j)+vL+u0+u(i,j+1))+-1*.25*dx*dy*-2*M*sin(M*i*dx)*cosh(M*j*dy);
                        
                    elseif(i==m) &&(j==n) %upper right corner check
                        u_new(i,j)= .25*(u_new(i-1,j)+vL+u_new(i,j-1)+uL)+-1*.25*dx*dy*-2*M*sin(M*i*dx)*cosh(M*j*dy);
                    
                    elseif (i==1) % x axis check
                        u_new(i,j)= .25*(v0+u(i+1,j)+u_new(i,j-1)+u(i,j+1))+-1*.25*dx*dy*-2*M*sin(M*i*dx)*cosh(M*j*dy);
                    
                    elseif (i==m)  % x=lx check
                        u_new(i,j)= .25*(u_new(i-1,j)+u_new(i,j-1)+u(i,j+1))+-1*.25*dx*dy*-2*M*sin(M*i*dx)*cosh(M*j*dy);
                    
                    elseif(j==1) % y axis check
                        u_new(i,j)= .25*(u_new(i-1,j)+u(i+1,j)+u(i,j+1))+-1*.25*dx*dy*-2*M*sin(M*i*dx)*cosh(M*j*dy);
                    
                    elseif(j==m) %y=ly
                        u_new(i,j)= .25*(u_new(i-1,j)+u(i+1,j)+u_new(i,j-1))+-1*.25*dx*dy*-2*M*sin(M*i*dx)*cosh(M*j*dy);
                    
                    else %all others wont depend on boundary
                        u_new(i,j)= .25*(u_new(i-1,j)+u(i+1,j)+u_new(i,j-1)+u(i,j+1))+-1*.25*dx*dy*-2*M*sin(M*i*dx)*cosh(M*j*dy);
                    
                    end
                     u(i,j)=u_new(i,j);
                end
           end
           u=u_new;
        end
        
     u
     solu
    

        
