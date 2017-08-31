                                       
clear;
clc;

% parameters
alpha = 0.30;    % capital share for CD tech.
beta  = 0.95;    % discount factor
delta = 0.10;    % depreciation rate

% construct K matrix: all possible k-combinations 
kmin=0.001;
kmax=10;
n=500;       % number of data points in the grid
ns= 3;        % number of values for the shock
kgrid = linspace(kmin,kmax,n)';     % column vector
K = kron(ones(n,1),kgrid');         % n number of row vectors
Kprime = K';                        % n number of column vectors

% space of shocks, and markov matrix 
A = [.5 1 1.5]; % Row vector !
Q = [.8 .1 .1; .1 .8 .1; .1 .1 .8];

% construct consumption matrix
tempc1=A(1)*(K.^alpha)+(1-delta)*K-Kprime;
tempc2=A(2)*(K.^alpha)+(1-delta)*K-Kprime;
tempc3=A(3)*(K.^alpha)+(1-delta)*K-Kprime;
tempc=[tempc1 tempc2 tempc3]; %nx3n
% remove invalid utility 
F=log(max(tempc,0)); %nx3n

%does the same thing with that in line 29
%tempc(tempc<0)=0;
%F=log(tempc);,,,,,,

iter= 0;
diff= 1;        % convergence criterion ------- renaming
tol = 0.0001;   % convergence parameter-------renaming   
v   = zeros(n,ns);

% value function iteration
while diff>tol;

    % construct the objective function
    Ev = v*Q; % expectation of value function, (nx3)(3x3)=nx3 matrix 
    vmatrix = [kron(ones(1,n),Ev(:,1)) kron(ones(1,n),Ev(:,2)) kron(ones(1,n),Ev(:,3))]; % nx3n
    
    [Tv,ig]=max(F+beta*vmatrix,[],1);    % max across columns; gives 1x3n row vector
    tempTv = reshape(Tv(1:3*n), n,[]); % nx3 matrix
    %tempTv=[Tv(1,1:n)' Tv(1,n+1:2*n)' Tv(1,2*n+1:3*n)'];   % nx3 matrix
    Tv=tempTv;

diff = max(max(abs(Tv-v))); % why do we need two max?
v = Tv;
display ([diff iter]) % -------------replacing
iter = iter+1;
end

ignew=reshape(ig(1:3*n),n,[]); %nx3 coloumn vector of ig
%tempkp=kron(ones(1,3),kgrid'); % 1x3n
%kp_row=tempkp(ig);
%kp=[kp_row(1,1:n)' kp_row(1,n+1:2*n)' kp_row(1,2*n+1:3*n)'];

kp = [kgrid(ignew(1:n,1)) kgrid(ignew(1:n,2)) kgrid(ignew(1:n,3))]; %nx3

cp=zeros(n,ns);
optc=zeros(n,1);
m=1;
i=1;
while i<=3,
   while m<=n,
   optc(m)= A(i)*kgrid(m)^alpha+(1-delta)*kgrid(m)-kp(m,i);
    m=m+1;
   end
   cp(:,i)=optc;
   i=i+1;
   m=1;
   end

%cp = [A(1)*kgrid.^alpha+(1-delta)*kgrid-kp(:,1) A(2)*kgrid.^alpha+(1-delta)*kgrid-kp(:,2) A(3)*kgrid.^alpha+(1-delta)*kgrid-kp(:,3)];

%Drawing optimal consumption policy, investment, next periods capital stock
%and value function
figure(1);
hold on;
subplot(2,2,1);
plot(kgrid, cp(:,1),'g', kgrid, cp(:,2),'b', kgrid, cp(:,3),'r');
hold off;
xlabel('current capital stock');
ylabel('optimal consumption');
legend('z_L','z_M','z_H','Location','southeast');


hold on;
subplot(2,2,2);
plot(kgrid, v(:,1),'g', kgrid, v(:,2),'b', kgrid, v(:,3),'r');
hold off;
xlabel('current capital stock');
ylabel('value function');
legend('z_L','z_M','z_H','Location','southeast');

hold on;
subplot(2,2,3);
plot(kgrid, kp(:,1),'g', kgrid, kp(:,2),'b', kgrid, kp(:,3),'r', kgrid, kgrid, 'k');
hold off;
xlabel('current capital stock');
ylabel('next periods capital stock');
legend('z_L','z_M','z_H','45 degree line','Location','southeast');

kcur=kron(ones(1,3), kgrid);
inv = kp - (1-delta)*kcur;

hold on;
subplot(2,2,4);
plot(kgrid, inv(:,1),'g', kgrid, inv(:,2),'b', kgrid, inv(:,3),'r');
hold off;
xlabel('current capital stock');
ylabel('optimal investment');
legend('z_L','z_M','z_H','Location','southeast');

z0 = A(1);
k0 = kmin;
T = 100;
p = rand (T,1);
z = zeros(T,1); % productivity shock vector
z(1) = z0;
for d = 2:100
    if p(d)<=0.8
        z(d)=z(d-1);
    elseif p(d)>.8 && p(d)<=0.9
            if z(d-1)~= A(1)
                z(d)=A(1);
            elseif z(d-1)~= A(2)
                z(d)=A(2);
            else 
                z(d)=A(3);
            end
    else
        if z(d-1)~= A(1)
                z(d)=A(1);
            elseif z(d-1)~= A(3)
                z(d)=A(3);
            else 
                z(d)=A(2);
        end
    end
end 

ct = zeros(T,1);
invt = zeros(T,1);
kpt = zeros(T,1);
kt = zeros(T,1);
yt = zeros(T,1);
pvt = zeros(T,1);
z_dummy = zeros(T,1);

for dd = 1:T
    if z(dd)==0.5
        z_dummy(dd)=1;
    elseif z(dd)==1
        z_dummy(dd)=2;
    else 
        z_dummy(dd)=3;
    end
end

kt(1) = kmin;
kpt(1) = kp(1,z_dummy(1));
yt(1) = z(1)*kt(1)^alpha;
ct(1) = yt(1)+(1-delta)*kt(1)-kpt(1);
invt(1)=kpt(1)-(1-delta)*kt(1);

for i = 2:T
kt(i)=kpt(i-1);
coloumn = z_dummy(i);
row = find(kgrid == kt(i));
kpt(i)=kp(row,coloumn);
yt(i)=z(i)*kt(i)^alpha;
ct(i)=yt(i)+(1-delta)*kt(i)-kpt(i);
invt(i)=kpt(i)-(1-delta)*kt(i);
pvt(i) = beta*v(row, coloumn);
end

%Drawing optimal consumption policy, investment, next periods capital stock
%and value function
time = linspace(1,100);

figure(2);

hold on;
subplot(2,3,1);
plot(time,z,'b');
hold off;
xlabel('time');
ylabel('sequence of shocks');
title('productivity shock') 

hold on;
subplot(2,3,2);
plot(time,kt,'b');
hold off;
xlabel('time');
ylabel('capital stock');
title('capital') 

hold on;
subplot(2,3,3);
plot(time,yt,'b');
hold off;
xlabel('time');
ylabel('output');
title('output') 

hold on;
subplot(2,3,4);
plot(time,ct,'g');
hold off;
xlabel('time');
ylabel('consumption');
title('consumption') 


hold on;
subplot(2,3,5);
plot(time,invt,'g');
hold off;
xlabel('time');
ylabel('investment');
title('investment') 

hold on;
subplot(2,3,6);
plot(time,pvt,'g');
hold off;
xlabel('time');
ylabel('PV of future expected discounted utility');
title('PV of future expected discounted utility') 

