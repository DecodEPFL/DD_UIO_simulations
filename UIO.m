clear

%% Parameters: 

n = 3;
m = 0;
md = 2;
p = 3;
Tini = 1;
%%% Tini MUST BE GREATER THAN OR EQUAL TO THE LAG OF THE SYSTEM!
Tf = 1;
MICROGRID=true;

%% System Matrices: 


if MICROGRID==true

    Rt = diag([2.0e-1, 3.0e-1, 1.0e-1, 5.0e-1]);
    Lt = diag([1.8e-3, 2.0e-3, 2.2e-3, 3.0e-3]);
    Ct = diag([2.2e-3, 1.9e-3, 1.7e-3, 2.5e-3]);

    K1 = diag([-2.134, -0.869, -0.480, -6.990]);
    K2 = diag([-0.163, -0.050, -0.108, -0.175]);
    K3 = diag([13.553, 48.285, 30.673, 102.960]);

    % Its = diag([4, 5, 6, 7]);

    alpha = (K1 - eye(4))/Lt;
    beta = (K2 - Rt)/Lt;
    gamma = K3/Lt;

    % delta = alpha;
    % 
    % R = diag([7.0e-2, 4.0e-2, 8.0e-2, 7.0e-2]);
    % L = diag([2.1e-6, 2.3e-6, 1.8e-6, 1.0e-6]);
    % 
    % YL = diag([1/9, 1/8, 1/10, 1/11]);
    
    IL = [6; 4; 3; 3.5];
    Vref = [47; 46.5; 48; 47.5];

    Ac = [        0         ,       1/Ct(1,1)    ,     0      ;
              alpha(1,1)    ,       beta(1,1)    , gamma(1,1) ;
                 -1         ,          0         ,     0      ];

    Bc = [];

    Ec = [ -1/Ct(1,1) , 0 ;
               0      , 0 ;
               0      , 1 ];

    Ts = 1e-2;

    A = expm(Ac*Ts);
    B = [];
    E = inv(Ac)*(A-eye(n))*Ec;

%     C = [1 0  zeros(1,n-2);
%          0 1  zeros(1,n-2)];
%     C = [0 1 0 zeros(1,n-3);
%          0 0 1 zeros(1,n-3)];
%     C = [1 0 0 zeros(1,n-3);
%          0 0 1 zeros(1,n-3)];
    C = eye(n);

else

    A = 2*rand(n,n)-1;
    A = A*rand/max(abs(eig(A)));
    % max(abs(eig(A)))

    B = 2*rand(n,m)-1;
    % B = [[0; 1; zeros(n-2,1)] [zeros(3,1); 1; zeros(n-4,1)]];

    E = 2*rand(n,md)-1;
    % E = [1; zeros(n-1,1)];

    C = 2*rand(p,n)-1;
    % C = [1  zeros(1,n-1)   ;
    %      0 0 1 zeros(1,n-3)];
    % C = [0 1  zeros(1,n-2)   ;
    %      0 0 0 1 zeros(1,n-4)];
    % C = [1  zeros(1,n-1)];
    % C = eye(n);

end

sys = ss(A, [B E], C, [], 1);

%% UIO Design

if rank(C*E) ~= rank(E)
    disp('UIO is not possible! Rank condition does not hold.')
%     return
end

H = E*pinv(E'*C')';
T = eye(size(H,1))-H*C;
% T*E

if rank(obsv(T*A,C)) ~= size(C,2)
    disp('UIO is not possible! Observability does not hold.')
%     return
end

%% Historical Data

Nhist = (Tini+Tf)*10;

x0 = 2*rand(n,1)-1;

ubar = 2*rand(m*Nhist,1)-1;
if MICROGRID==true
    dbar = 1*(2*rand(md*Nhist,1)-1)+repmat([IL(1)+5; Vref(1)+3],Nhist,1);
else
    dbar = 1*(2*rand(md*Nhist,1)-1);
end

[ysim,time,xsim] = lsim(sys,[reshape(ubar,m,Nhist)' reshape(dbar,md,Nhist)'],0:Nhist-1,x0);

ybarvec = reshape(ysim',Nhist*p,1);
xbarvec = reshape(xsim',Nhist*n,1);

% D = compute_Hankel(dbar, md, Tini+Tf);
% Dp = D(1:Tini*md,:);
% Df = D(Tini*md+1:end,:);

Y = compute_Hankel(ybarvec, p, Tini+Tf);
Yp = Y(1:Tini*p,:);
Yf = Y(Tini*p+1:end,:);

X = compute_Hankel(xbarvec, n, Tini+Tf);
Xp = X(1:Tini*n,:);
Xf = X(Tini*n+1:end,:);

if m~=0
    U = compute_Hankel(ubar, m, Tini+Tf);
    Up = U(1:Tini*m,:);
    Uf = U(Tini*m+1:end,:);
else
    U = [];
    Up = [];
    Uf = [];
end
% %% Recent Data
% 
% x0rec = 2*rand(n,1)-1;
% 
% usim = 2*rand(m*(Tini+Tf),1)-1;
% dsim = 10*(2*rand(md*(Tini+Tf),1)-1);
% 
% [ysim,time,xsim] = lsim(sys,[reshape(usim,m,(Tini+Tf))' reshape(dsim,md,(Tini+Tf))'],0:(Tini+Tf)-1,x0rec);
% 
% yvec = reshape(ysim',(Tini+Tf)*p,1);
% xvec = reshape(xsim',(Tini+Tf)*n,1);
% 
% uini = usim(1:Tini*m);
% yini = yvec(1:Tini*p);
% xini = xvec(1:Tini*n);
% % xini = zeros(Tini*n,1);
% 
% u = usim(Tini*m+1:end);
% atk = rand(Tf*p,1);
% atk = zeros(Tf*p,1);
% y = yvec(Tini*p+1:end)+atk;
% x = xvec(Tini*n+1:end);
% 
% xpred = Xf*pinv([Up;Yp;Xp;Uf;Yf])*[uini;yini;xini;u;y];
% 
% % plot(x)
% % hold on
% % plot(xpred)
% 
% % max(abs(x-xpred))
% max(abs(y-kron(eye(Tf),C)*xpred))
% 
% % figure()
% % plot(x(1:n:end))
% % hold on 
% % plot(xpred(1:n:end))


%% Online Prediction (one-step prediction at each point)

x0rec = [46.78; 9; 11]+1*(2*rand(n,1)-1);

Tsim = (Tini+Tf)*50;  % Length of the simulation

usim = 2*rand(m*Tsim,1)-1;
if MICROGRID==true
    dsim = .1*(2*rand(md*Tsim,1)-1)+repmat([IL(1)+3; Vref(1)],Tsim,1);
else
    dsim = 1*(2*rand(md*Tsim,1)-1);
end

[ysim,time,xsim] = lsim(sys,[reshape(usim,m,Tsim)' reshape(dsim,md,Tsim)'],0:Tsim-1,x0rec);

yvec = reshape(ysim',Tsim*p,1);
xvec = reshape(xsim',Tsim*n,1);

uini = usim(1:Tini*m);
yini = yvec(1:Tini*p);
xini = xvec(1:Tini*n);
% xini = zeros(Tini*n,1);

uf = usim(Tini*m+1:end);
atkidx = 50*Tini;
atk = zeros((Tsim-Tini)*p,1);
% atk((atkidx-Tini)*p+1:end) = .1*(2*rand((Tsim-atkidx)*p,1)-1);
atk((atkidx-Tini)*p+1:end) = .1*(ones((Tsim-atkidx)*p,1));
yf = yvec(Tini*p+1:end)+atk;
xf = xvec(Tini*n+1:end);

u = [uini; uf];
y = [yini; yf];
x = [xini; xf];

xpred = 10*rand(Tsim*n,1);
% xpred(1:n) = [42; 5; 0]+1*(2*rand(n,1)-1);
xpred(1:n) = x0rec+.1*(2*rand(n,1)-1);
W = pinv([Up;Yp;Xp;Uf;Yf]);
W1 = W(:, 1:(p+m)*Tini);
W2 = W(:, (p+m)*Tini+1:(p+m+n)*Tini);
W3 = W(:, (p+m+n)*Tini+1:end);
abs(eig(Xf*W2))

if max(abs(eig(Xf*W2)))>=1
    disp('Data-driven UIO is not possible! Stability condition does not hold.')
end
if rank([null([Up;Yp;Xp;Uf;Yf]), null(Xf)]) ~= rank(null(Xf))
    disp('Data-driven UIO is not possible! Null space condition does not hold.')
end

for t = Tini:(Tsim-1)
    xpredini = xpred((t-Tini)*n+1:t*n);
    uoini = u((t-Tini)*m+1:t*m);
    yoini = y((t-Tini)*p+1:t*p);
    uof = u(t*m+1:t*m+m);
    yof = y(t*p+1:t*p+p);
    xpred(t*n+1:t*n+Tf*n) = Xf*W*[uoini;yoini;xpredini;uof;yof];
%     xpred(t*n+1:t*n+Tf*n) = Xf*lsqminnorm([Up;Yp;Xp;Uf;Yf], [uini;yini;xpredini;u;y]);
end

ypred = kron(eye(Tsim),C)*xpred;

if MICROGRID==true
    figure()
    subplot(3,1,1)
    plot(time, x(1:n:end), '-b', 'linewidth', 2)
    hold on 
    grid on
    plot(time, xpred(1:n:end), '--r', 'linewidth', 2)
    yl = ylim;
    plot([atkidx atkidx],ylim,':k', 'linewidth', 2)
    ylabel('Voltage ($V$)', 'interpreter', 'latex', 'fontsize', 16)
    legend({'Real','Estimate'}, 'interpreter', 'latex', 'fontsize', 16)
    subplot(3,1,2)
    plot(time, x(2:n:end), '-b', 'linewidth', 2)
    hold on 
    grid on
    plot(time, xpred(2:n:end), '--r', 'linewidth', 2)
    yl = ylim;
    plot([atkidx atkidx],ylim,':k', 'linewidth', 2)
    ylabel('Current ($A$)', 'interpreter', 'latex', 'fontsize', 16)
    subplot(3,1,3)
    plot(time, x(3:n:end), '-b', 'linewidth', 2)
    hold on 
    grid on
    plot(time, xpred(3:n:end), '--r', 'linewidth', 2)
    yl = ylim;
    plot([atkidx atkidx],ylim,':k', 'linewidth', 2)
    ylabel('Integrator ($Vs$)', 'interpreter', 'latex', 'fontsize', 16)
    xlabel('Timestep $t$', 'interpreter', 'latex', 'fontsize', 16)
    
    figure()
    subplot(3,1,1)
    plot(time, (x(1:n:end)-xpred(1:n:end)), '-b', 'linewidth', 2)
    hold on
    grid on
    yl = ylim;
    plot([atkidx atkidx],ylim,':k', 'linewidth', 2)
    ylabel('$e_{t,1}$ ($V$)', 'interpreter', 'latex', 'fontsize', 16)
%     ylabel('$V_i-\hat{V}_i$ ($V$)', 'interpreter', 'latex', 'fontsize', 16)
    subplot(3,1,2)
    plot(time, (x(2:n:end)-xpred(2:n:end)), '-b', 'linewidth', 2)
    hold on
    grid on
    yl = ylim;
    plot([atkidx atkidx],ylim,':k', 'linewidth', 2)
    ylabel('$e_{t,2}$ ($A$)', 'interpreter', 'latex', 'fontsize', 16)
%     ylabel('$I_{ti}-\hat{I}_{ti}$ ($A$)', 'interpreter', 'latex', 'fontsize', 16)
    subplot(3,1,3)
    plot(time, (x(3:n:end)-xpred(3:n:end)), '-b', 'linewidth', 2)
    hold on
    grid on
    yl = ylim;
    plot([atkidx atkidx],ylim,':k', 'linewidth', 2)
    ylabel('$e_{t,3}$ ($Vs$)', 'interpreter', 'latex', 'fontsize', 16)
%     ylabel('$v_i-\hat{v}_i$ ($Vs$)', 'interpreter', 'latex', 'fontsize', 16)
    xlabel('Timestep $t$', 'interpreter', 'latex', 'fontsize', 16)
    
    figure()
    subplot(3,1,1)
    plot(time, (y(1:p:end)-ypred(1:p:end)), '-b', 'linewidth', 2)
    hold on
    grid on
    yl = ylim;
    plot([atkidx atkidx],ylim,':k', 'linewidth', 2)
    ylabel('$r_{t,1}$ ($V$)', 'interpreter', 'latex', 'fontsize', 16)
%     ylabel('$V_i-\hat{V}_i$ ($V$)', 'interpreter', 'latex', 'fontsize', 16)
    subplot(3,1,2)
    plot(time, (y(2:p:end)-ypred(2:p:end)), '-b', 'linewidth', 2)
    hold on
    grid on
    yl = ylim;
    plot([atkidx atkidx],ylim,':k', 'linewidth', 2)
    ylabel('$r_{t,2}$ ($A$)', 'interpreter', 'latex', 'fontsize', 16)
%     ylabel('$I_{ti}-\hat{I}_{ti}$ ($A$)', 'interpreter', 'latex', 'fontsize', 16)
    subplot(3,1,3)
    plot(time, (y(3:p:end)-ypred(3:p:end)), '-b', 'linewidth', 2)
    hold on
    grid on
    yl = ylim;
    plot([atkidx atkidx],ylim,':k', 'linewidth', 2)
    ylabel('$r_{t,3}$ ($Vs$)', 'interpreter', 'latex', 'fontsize', 16)
%     ylabel('$I_{ti}-\hat{I}_{ti}$ ($A$)', 'interpreter', 'latex', 'fontsize', 16)
    xlabel('Timestep $t$', 'interpreter', 'latex', 'fontsize', 16)
else
    figure()
    subplot(2,2,1)
    plot(time, x(1:n:end), '-b')
    hold on 
    plot(time, xpred(1:n:end), '--r')
    ylabel('x1')
    subplot(2,2,2)
    plot(time, x(2:n:end), '-b')
    hold on 
    plot(time, xpred(2:n:end), '--r')
    ylabel('x2')
    subplot(2,2,3)
    plot(time, x(3:n:end), '-b')
    hold on 
    plot(time, xpred(3:n:end), '--r')
    ylabel('x3')
    subplot(2,2,4)
    plot(time, x(4:n:end), '-b')
    hold on 
    plot(time, xpred(4:n:end), '--r')
    ylabel('x4')

    figure()
    subplot(2,2,1)
    plot(time, (x(1:n:end)-xpred(1:n:end)), '-b')
    ylabel('e1')
    subplot(2,2,2)
    plot(time, (x(2:n:end)-xpred(2:n:end)), '-b')
    ylabel('e2')
    subplot(2,2,3)
    plot(time, (x(3:n:end)-xpred(3:n:end)), '-b')
    ylabel('e3')
    subplot(2,2,4)
    plot(time, (x(4:n:end)-xpred(4:n:end)), '-b')
    ylabel('e4')

end

%% FUNCTION TO COMPUTE HANKEL MATRIX FROM DATA SEQUENCE
function H = compute_Hankel(w,dim,L)
    Nd = floor(size(w,1)/dim);
    H = [];
    for i = 1:Nd-L+1
        H(:,i) = w((i-1)*dim+1:(i+L-1)*dim);
    end
end

