clc,clear
T = 1 ;
S_d = 30 ;
S_u = 83 ;
X_d = log(S_d) ;
X_u = log(S_u) ;
K = 50 ;
r = 0.08 ;
sigma = 0.25 ;
alpha = 1.5 ;

N = 200 ;
M = 300 ;
dt = T/N ;
dx = (X_u-X_d)/M ;
t = linspace(0,T,N+1) ;
X = linspace(X_d,X_u,M+1) ;
S = exp(X) ;

v = -sigma^alpha/cos(alpha*pi/2)/2 ;
zeta = v*dt/(2*dx^alpha) ;
xi = dt*(r-v)/(4*dx) ;
eta = dt*r/2 ;

g = zeros(M+1,1);
omega = zeros(M+1,1);
g(1,1) = 1 ;
omega(1,1) = alpha*g(1,1)/2 ;
for k = 2:M+1
    g(k,1) = (1-(alpha+1)/(k-1))*g(k-1,1);
    omega(k,1) = alpha*g(k,1)/2+(2-alpha)*g(k-1,1)/2 ;
end
A = zeros(M-1,M-1);
for i = 1:M-1
    for j = 1:M-1
        if i==j
            A(i,j) = omega(2,1);
        elseif i-j == 1
            A(i,j) = omega(3,1);
        elseif i-j == -1
            A(i,j) = omega(1,1);
        elseif i-j >= 2
            A(i,j) = omega(i-j+2,1);
        else
            A(i,j) = 0;
        end
    end
end
B = diag(-1*ones(M-2,1),-1)+diag(ones(M-2,1),1) ;
L = (1+eta)*eye(M-1)-zeta*A-xi*B ;
R = (1-eta)*eye(M-1)+zeta*A+xi*B ;

V = zeros(N+1,M+1);
for i = 1:M+1
    if X_d+(i-1)*dx <= log(83)
        V(1,i) = max(exp(X_d+(i-1)*dx)-K,0);
    else
        V(1,i) = 0;
    end
end

V(:,1) = zeros(N+1,1);
V(:,M+1) = zeros(N+1,1);
for i = 2:N+1
    V(i,2:M) = (L\R*(V(i-1,2:M))')';
end

%surf(S,T-t,V,'EdgeColor','none') 
%xlabel('Price of stock'); 
%ylabel('Time') 
%zlabel('Price of option') 

plot(S,V(1,:));hold on;
plot(S,V(10,:));hold on;
plot(S,V(20,:));hold on;
plot(S,V(30,:));hold on;
plot(S,V(40,:));hold on;
plot(S,V(50,:));
xlabel('Price of stock')
ylabel('Price of option')
handle = legend('$t=0$','$t=\frac{1}{20}$','$t=\frac{1}{10}$','$t=\frac{3}{20}$','$t=\frac{1}{5}$','$t=\frac{1}{4}$','location','northwest');
set(handle,'Interpreter','latex', 'FontSize',12)


