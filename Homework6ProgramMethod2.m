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
%MAKE SURE THAT IT IS EVEN AND DIVISIBLE BY 4; use 20 as minimum
m=160;
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

iteration=1000;
%defining x & y intervals

x=[0:dx:lx];
y=[0:dy:ly];

%Closed-form solution given in Assignment sheet 
% M an integer
M=6;
for p=1:length(x)
    for q=1:length(y)
      solU(p,q)=(ly-y(q)).*sin(M.*x(p)).*sinh(M.*y(q));
      % given f in del2u=f
      solF(p,q)= -2.*M.*sin(M.*x(p)).*cosh(M.*y(q));
    end
end




%setting up u properly: i.e. adding the boundary conditions
u(m+2,n+2)=0; % will bc's on the 'plate'
u(:,1)=u0;
u(:,m+2)=uL;
u(1,:)=v0;
u(n+2,:)=vL;

L1_error=0;
error=0;
u_new(m+2,n+2)=0;

%This nested loop solves for the first u values
% and then applies the gauss-seidel to come w/
%better values at higher iterations
%to not compute too much in a loop, will do the following 
beta= -dx*dy*2*M;
for s=1:iteration
    for j=2:n+1 %done this way to ignore the bc entries on the matrix
        for i=2:m+1 % see above
       
       u(i,j)= (beta*sin(M*(i-1)*dx)*cosh(M*(j-1)*dy) -u_new(i,j-1)-u_new(i-1,j)-u(i+1,j) -u(i,j+1))/-4;
        u_new(i,j)=u(i,j); % new value assigned as in gauss-seidel
        end %end of i loop
        
    end % end of j loop
   u=u_new;      %iteration happens here
end

%Calculating the L1 error
     for i=2:m+1
         for j=2:n+1
             error= abs(u(i,j)-solU(i,j));
           L1_error=L1_error+error;
         end
     end
     L1_error= L1_error/(m*n)
% %Here prints out the surface for the numerical and closed form
% %Uncomment to see when program is compiles
%    % be sure to comment 2-d plots below
   k=surf(x,y,solU); %solU is closed-form solution 
     %makes true solution a little transparent so can see both surfaces
     %easier
    alpha(k,.4)
   hold on
   
   surf(x,y,u)
   title(['Numerical vs Closed-form surfaces when the m is= ' num2str(m)])
   
                    %plotting on 2-d points
                    %comment everything below is doing surface plot
    %setting up points to judge closeness in 2-d at
    % can choose other points but be sure that all are whole numbers
length=length(solU);
quarter= (length-2)*.5 +1 
crescent= (length-2)*.25 +1
gibbous= (length-2)*.75+1
       
                %AT Xs
uTruehalfx(length)=0;
uTruehalfx=solU(quarter,:);
uNumhalfx(length)=0;
uNumhalfx=u(quarter,:);
diff1=uNumhalfx-uTruehalfx;

uTruegibbousx(length)=0;
uTruegibbousx=solU(gibbous,:);
uNumgibbousx(length)=0;
uNumgibbousx=u(gibbous,:);
diff3= uNumgibbousx-uTruegibbousx;

uTruequarterx(length)=0;
uNumquarterx(length)=0;
uTruequarterx=solU(crescent,:);
uNumquarterx=u(crescent,:);
diff4=uNumquarterx-uTruequarterx;

figure

hold on
plot(x,uNumhalfx,'+g')
plot(x,uTruehalfx,'g')
plot(x,uNumquarterx,'+r')
plot(x,uTruequarterx,'r')
plot(x,uTruegibbousx,'c')
plot(x,uNumgibbousx,'+c')
legend('numerical half-way', 'closed-form half-way', 'numerical 1/4 ', 'closed-form at 1/4','numerical 3/4','closed-form 3/4')
title([' Closed-form and numercial values at various xs: 3 xs selected and "y" is allowed to vary  on them. iterations= ' num2str(iteration) '  mesh size=' num2str(m)])
xlabel(' y values')
ylabel('u values')
% At Y's
uTruehalfy(length)=0;
uTruehalfy=solU(:,quarter);
uNumhalfy(length)=0;
uNumhalfy=u(:,quarter);
diff1y=uNumhalfy-uTruehalfy;

uTruegibbousy(length)=0;
uTruegibbousy=solU(:,gibbous);
uNumgibbousy(length)=0;
uNumgibbousy=u(:,gibbous);
diff3y= uNumgibbousy-uTruegibbousy;

uTruequartery(length)=0;
uNumquartery(length)=0;
uTruequartery=solU(:,crescent);
uNumquartery=u(:,crescent);
diff4y=uNumquartery-uTruequartery;

figure

hold on
plot(y,uNumhalfy,'+g')
plot(y,uTruehalfy,'g')
plot(y,uNumquartery,'+r')
plot(y,uTruequartery,'r')
plot(y,uTruegibbousy,'c')
plot(y,uNumgibbousy,'+c')
legend('numerical half-way', 'closed-form half-way', 'numerical 1/4 ', 'closed-form at 1/4','numerical 3/4','closed-form 3/4')
title([' Closed-form and numercial values at various ys: 3 ys selected and "x" is allowed to vary  on them. iterations= ' num2str(iteration) '  mesh size=' num2str(m)])

xlabel(' x values')
ylabel('u values')