function Example1

[nx, nu, nw, ny, nz, r, A, B, E, C, D, G, H, J] = SysParas_Ex1;

T  = 50;
Tsamp = 2;
N  = T/Tsamp + 1;

x  = zeros(nx,N); x(:,1)  = [0.4;-0.3];
u  = zeros(nu,N);
xh = zeros(nx,N); xh(:,1) = [0; 0];
y  = zeros(ny,N);
z  = zeros(nz,N);
th1  = zeros(1,N);
th1h = zeros(1,N);
th2  = zeros(1,N);
th2h = zeros(1,N);



Hinf = zeros(1,N);
Zsum = zeros(1,N);
Wsum = zeros(1,N);


[gama, L_, F_, S, Z] = LMIs_Ex1(0.2,0.4);


for k = 1:N-1
    %% Disturbance
    wk = exp(-0.5*k)*sin(0.5*k);

    %% Fuzzy basic functions in Plant Side
    th1(k) = (1 - sin(x(1,k)))/2;
    th2(k) = 1 - th1(k);

    %% Fuzzy basic functions in Controller Side
    th1h(k) = (1 - sin(xh(1,k)))/2;
    th2h(k) = 1 - th1(k);
      
    
    %% System matrices and output
    Ath  =  th1(k)*A(:,:,1)  +  th2(k)*A(:,:,2);
    Bth  =  th1(k)*B(:,:,1)  +  th2(k)*B(:,:,2);
    Eth  =  th1(k)*E(:,:,1)  +  th2(k)*E(:,:,2);
    Cth  =  th1(k)*C(:,:,1)  +  th2(k)*C(:,:,2);
    Dth  =  th1(k)*D(:,:,1)  +  th2(k)*D(:,:,2);
    Jth  =  th1(k)*J(:,:,1)  +  th2(k)*J(:,:,2);
    Gth  =  th1(k)*G(:,:,1)  +  th2(k)*G(:,:,2);
    Hth  =  th1(k)*H(:,:,1)  +  th2(k)*H(:,:,2);
    
      
    %% Feedback and Observer Gains and Matrices
    Athh  =  th1h(k)*A(:,:,1)  +  th2h(k)*A(:,:,2);
    Bthh  =  th1h(k)*B(:,:,1)  +  th2h(k)*B(:,:,2);
    Cthh  =  th1h(k)*C(:,:,1)  +  th2h(k)*C(:,:,2);
    
    
    Fthh  =  inv(th1h(k)*S(:,:,1) + th2h(k)*S(:,:,2))*(th1h(k)*F_(:,:,1) + th2h(k)*F_(:,:,2));
    Lthh  =  inv(th1h(k)*Z(:,:,1) + th2h(k)*Z(:,:,2))*(th1h(k)*L_(:,:,1) + th2h(k)*L_(:,:,2));
    
   
    
    u(k)   =  Fthh*xh(:,k);
    %% Observer
    y(:,k)   = Cth*x(:,k) + Dth*wk;
    xh(:,k+1) = Athh*xh(:,k) + Bthh*u(k) + Lthh*(y(:,k) - Cthh*xh(:,k));
    
    %% System
    x(:,k+1)  = Ath*x(:,k) + Bth*u(k) + Eth*wk;
    z(k)      = Gth*x(:,k) + Hth*u(k) + Jth*wk;
    
    Zsum(k+1) = Zsum(k) + z(k)'*z(k);
    Wsum(k+1) = Wsum(k) + wk'*wk;
    Hinf(k+1) = sqrt(Zsum(k+1)/(gama^2*Wsum(k+1)));
end

fontsize = 14;
linewidth = 2;

clf(figure(11)); axes('Position',[0.1 0.1 0.87 0.87]);
stairs(0:N-1,x(1,:),'linewidth',linewidth); hold on;
stairs(0:N-1,xh(1,:),'linewidth',linewidth);
set(gca,'fontsize',fontsize);
axis([0 N-1 -0.1 0.45]);
legend('$x_{1,k}$','$\hat x_{1,k}$','fontsize',fontsize+6,'interpreter','latex')
grid on;

clf(figure(12));  axes('Position',[0.1 0.1 0.87 0.87]);
stairs(0:N-1,x(2,:),'linewidth',linewidth); hold on;
stairs(0:N-1,xh(2,:),'linewidth',linewidth);
set(gca,'fontsize',fontsize);
axis([0 N-1 -0.5 0.45]);
legend('$x_{2,k}$','$\hat x_{2,k}$','fontsize',fontsize+6,'interpreter','latex')
grid on;

clf(figure(14)); axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,u,'linewidth',linewidth);
legend('$u(k)$','fontsize',fontsize,'interpreter','latex');
set(gca,'fontsize',fontsize);
% axis([0 N-1 -0.8 1.1]);
grid on;

% clf(figure(7));
% stairs(0:N-1,th1-th1h,'linewidth',linewidth);
% set(gca,'fontsize',fontsize);
% legend('$|\theta_i(k) - \hat\theta_i(k)|$','fontsize',fontsize,'interpreter','latex');
% 
% clf(figure(5)); axes('Position',[0.1 0.1 0.85 0.85]);
% stairs(0:N-1,Hinf,'linewidth',linewidth);
% set(gca,'fontsize',fontsize);
% axis([0 N-1 0 1]);
% legend('$\frac{\sum_{k=0}^{T}\Vert z_k \Vert_2^2}{\gamma_{min}^2 \sum_{k=0}^{T}\Vert w_k \Vert_2^2}$'...
%     ,'fontsize',fontsize+2,'interpreter','latex','location','best');
% grid on;

end




