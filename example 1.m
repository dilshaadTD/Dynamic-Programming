%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Macro Pset #4 - Q4 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

clear;
clc;

% parameters
sigma = 1.50;    
beta  = 0.99322;    % discount factor
w = 1;   %income when e=1
b = 0.1; %unemployment benefit when e=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Bisection Method to find stationary interest rate r;
% In complete market, beta*(1+r)=1, so proper guess r>0.00683;
r= -0.0105;        %initial guess;
%rmax=0.1; 
%rmin=-0.1;
tol_r=0.005; %given
diff_r=1;
iter_r=0;
while (diff_r>tol_r),
    
%continue w/ given guess r until checking mrk clr, & then update if needed;
   % r=(rmax+rmin)/2;
   
   
% construct agrid for with shock e=0,1;
amin= -2;
amax=2;     %adjust if needed
n=300;       %number of grid pts
ns= 2;       %number of realizations for shock e=0,1
agrid = linspace(amin,amax,n)';     % nx1
A = kron(ones(n,1),agrid');         % n number of row vectors of agrid
Aprime = A';                        % n number of column vectors of agrid

e = [1 0]; % Idiosyncratic shock on employment status
Q = [0.9250 0.5; 0.075 0.5]; %Transition Matrix for e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% solve HH problem by Value Function Iteration 
% first, construct consumption matrix
tempc1 = (1+r)*A + e(1)*w + (1-e(1))*b - Aprime; %tempc1 when e=1;
tempc2 = (1+r)*A + e(2)*w + (1-e(2))*b - Aprime; %tempc1 when e=0;
%c1=tempc1.^(1-sigma)/(sigma).*(tempc1>=0)-10^10*(tempc1<0);
%c2=tempc2.^(1-sigma)/(sigma).*(tempc2>=0)-10^10*(tempc2<0);
tempc=[tempc1 tempc2]; %nx2n 

F=((max(tempc,0)).^(1-sigma))/(1-sigma); %nx2n

% start value function iteration;
iter= 0;
diff= 1;       
tol = 0.0001;  
v   = zeros(n,ns); %nx2: initial guess for v

while diff>tol;

    % construct the objective function
    Ev = v*Q; % expectation of value function, (nx2)(2x2)=nx2 matrix 
    vmatrix = [kron(ones(1,n),Ev(:,1)) kron(ones(1,n),Ev(:,2))]; % nx2n
    
    [Tv,ig]=max(F+beta*vmatrix,[],1);  % 1x2n: max for each column; 
    tempTv = reshape(Tv(1:2*n), n,[]); % nx2: shock e=1 for 1st column;
    Tv=tempTv;

diff = max(max(abs(Tv-v)));
v = Tv;
%display ([diff iter])
iter = iter+1;
end

% policy function for a'
ignew=reshape(ig(1:2*n),n,[]); %1x2n => nx2; 
ap = [agrid(ignew(1:n,1)) agrid(ignew(1:n,2))]; %nx2
kl;
% optimal consumption matrix under ap;
cp=zeros(n,ns);
optc=zeros(n,1);
m=1;
i=1;
while i<=2,
   while m<=n,
   optc(m)= max((1+r)*agrid(m)+e(i)*w+(1-e(i))*b-ap(m,i),0);
    m=m+1;
   end
   cp(:,i)=optc;
   i=i+1;
   m=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% constructing H mapping for lambda using Q and ap;
% fisrt, construct indicator function for a'=a'(e,a);
indicator=zeros(n,n,2);
%B=zeros(2*n,2*n);

for j=1:2
    for k=1:n
        indicator(k,ignew(k,j),j)=1;
    end
    B(j+(j-1)*(n-1):j*n,:)=kron(Q(j,:),indicator(:,:,j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% calculate stationary distribution, lambda_star by iteration 
tol_l=0.0001;
diff_l=1;
iter_l=0;
lambda=kron(ones(2*n,1),1/(2*n));   % 2nx1
B=B'; %now today's state is given across columns;

while(diff_l>tol_l) %&& (iter_l<=500),
    
    Tlambda=B*lambda;  %(2nx2n)(2nx1)  
    diff_l=max(abs(Tlambda-lambda));
    lambda=Tlambda; %initial guess for lambda_star 
    iter_l=iter_l+1;
    %display([diff_l iter_l])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%check market clearing condition;
%Tlambdanew=reshape(Tlambda,n,[])'; %(2nx1)=>(nx2)=>(2xn)

apnew=reshape(ap,1,[]);
%excess_a=Tlambda(1:n)'*agrid+Tlambda(n+1:2*n)'*agrid;
excess_a=0;
i = 0;
for i=1:2*n;
    if i<=n
        excess_a=excess_a+Tlambda(i)*agrid(i);
    elseif i>n
        excess_a=excess_a+Tlambda(i)*agrid(i-n);
    end
end
        
%Tlambda*apnew; %(2nx1)(1x2n)
if excess_a>0
   % rmax=r; % if excess supply, decrease r;
r=r-0.0001;
else
  %  rmin=r; % if excess demand, increase r;
r=r+0.0001;
end

diff_r=abs(excess_a);
iter_r=iter_r+1;

disp([excess_a r iter_r]);
end

yrate=(r+1)^6-1;    % annual interest rate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% plot value function v(a,e), contingent on the e=0,1;
subplot(2,2,1)
plot(agrid,Tv(1:n,1),'r')
hold on
plot(agrid,Tv(1:n,2),'b')
title('value functions')
xlabel('a')
ylabel('v(a,e)')
legend('employed (e=1)','unemployed (e=0)','Location', 'southeast')

% plot policy function for consumption c(a,e), contingent on the e=0,1;
subplot(2,2,2)
plot(agrid,cp(1:n,1),'r')
hold on
plot(agrid,cp(1:n,2),'b')
title('policy functions for consumption c(a,e)')
xlabel('a')
ylabel('c(a,e)')
legend('employed (e=1)','unemployed (e=0)','Location', 'northwest')

% plot next period asset level a'(a,e), contingent on the e=0,1;
subplot(2,2,3)
plot(agrid,ap(1:n,1),'r')
hold on
plot(agrid,ap(1:n,2),'b')
hold on
plot(agrid,agrid,'g')
title('policy functions for next period asset')
xlabel('a')
ylabel('policy function for aprime')
legend('employed (e=1)','unemployed (e=0)','45 Degree Line','Location', 'northwest')

% plot stationary distribution lambda(a,e), contingent on the e=0,1;
subplot(2,2,4)
plot(agrid,Tlambda(1:n,1),'r')
hold on
plot(agrid,Tlambda(n+1:2*n,1),'b')
title('statrionary distribution')
xlabel('a')
ylabel('lambda(a,e) (distribution)')
legend('employed (e=1)','unemployed (e=0)','Location', 'northwest')

