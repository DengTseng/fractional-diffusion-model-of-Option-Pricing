alpha = [1.3,1.5,1.7,2];
S_d = 0.1 ;
S_u = 100 ;
X_d = log(S_d) ;
X_u = log(S_u) ;
M = 300 ;
X = linspace(X_d,X_u,M+1) ;
S = exp(X) ;
for a = alpha
    plot(S,fractional2(a)-fractional2(2));hold on;
end
xlabel('Price of stock');
handle = ylabel('$V_{FMLS}-V_{BS}$');
set(handle,'Interpreter','latex')
handle2 = legend('$\alpha = 1.3$','$\alpha = 1.5$','$\alpha = 1.7$','$\alpha = 2.0$','location','northwest');
set(handle2,'Interpreter','latex', 'FontSize',12)