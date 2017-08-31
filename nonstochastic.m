clear;
clc;

%parameters
global beta sigma 
alpha=0.33;
beta=.96;
delta=0.1;

%steady state
ks=(alpha/(1/beta-(1-delta)))^(1/(1-alpha));
cs = ks^alpha -delta*ks;


kl=0.85*ks;
kh=1.15*ks;
n=150;
k=linspace(kl,kh,n);


k=k'; % k vectorunu colomn cinsinden yazýyorum
K=kron(ones(1,n),k);
c=K'.^alpha+(1-delta)*K'-K;%consumptions


%sigma=1;

C=zeros(n,4);
KK=zeros(n,4);

for sigma = 1:4 
    U=utilfun(c);

% Value function iteration
tol=10^-8;
diff=1;
Tv=ones(n,1);
   while diff>tol;
     v=Tv;
     vmat=kron(ones(1,n),v);
     [Tv,ig]=max(U+beta*vmat);
     Tv=Tv';
     diff=max(abs(v-Tv));
     disp (diff);
    end

kp = k(ig); %policy function
optc=zeros(n,1); %optimal consýmption decision
 
m=1;
   while m<n,
   optc(m+1)= kp(m)^alpha+(1-delta)*kp(m)-kp(m+1);
   %optc(m)=k(m)^alpha + (1-delta)*k(m) - kp(m);
    m=m+1;
   end

C(:,sigma)=optc;
KK(:,sigma)=kp;
end

hold on;
subplot(1,2,1);
plot(k,KK(:,1),'g-',k,KK(:,2),'r--',k,KK(:,3),'y:',k,KK(:,4),'b-.',k,k,':');
hold off;


xlabel('current capital stock');
ylabel('next period capital stock');
title('g(k(t))=k(t+1)');
legend('sigma=1','sigma=2','sigma=3','sigma=4','45 degree line', 'Location','northwest');

hold on;
subplot(1,2,2);
plot(k,C(:,1),'g-',k,C(:,2),'r--',k,C(:,3),'y:',k,C(:,4),'-.');
hold off;

xlabel('current capital stock');
ylabel('current consumption');
title('optimal consumption decision');
legend('sigma=1','sigma=2','sigma=3','sigma=4','Location','northeast');


    
    

