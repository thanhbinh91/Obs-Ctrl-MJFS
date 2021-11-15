function Example2

[nx, nu, nw, ny, nz, s, r, A, B, E, C, D, G, H, J, Hk, Hb, Up_Pi, Lo_Pi] = SysParas_Ex2;

T  = 50;
Tsamp = 2;
N  = T/Tsamp + 1;

x  = zeros(nx,N); x(:,1)  = [0.5; -0.2];
u  = zeros(nu,N);
xh = zeros(nx,N); xh(:,1) = [0; 0];
y  = zeros(ny,N);
z  = zeros(nz,N);
th1  = zeros(1,N);
th1h = zeros(1,N);
th2  = zeros(1,N);
th2h = zeros(1,N);
rho  = ones(1,N);


Hinf = zeros(1,N);
Zsum = zeros(1,N);
Wsum = zeros(1,N);


[gama, L_, F_, S, Z] = LMIs_Ex2(0.1,0.4);

Pi    = [.4  .4  .2;
         .3  .6  .1;
         .5  .2  .3];

for k = 1:N-1
    wk = exp(-0.5*k)*sin(0.5*k);
    g = rho(k);

    %% Fuzzy basic functinos in Plant Side
    th1(k) = (1 - sin(x(1,k)))/2;
    th2(k) = 1 - th1(k);

    %% Fuzzy basic functinos in Controller Side
    th1h(k) = (1 - sin(xh(1,k)))/2;
    th2h(k) = 1 - th1(k);
      
    
    %% System matrices and output
    Ath  =  th1(k)*A(:,:,g,1)  +  th2(k)*A(:,:,g,2);
    Bth  =  th1(k)*B(:,:,g,1)  +  th2(k)*B(:,:,g,2);
    Eth  =  th1(k)*E(:,:,g,1)  +  th2(k)*E(:,:,g,2);
    Cth  =  th1(k)*C(:,:,g,1)  +  th2(k)*C(:,:,g,2);
    Dth  =  th1(k)*D(:,:,g,1)  +  th2(k)*D(:,:,g,2);
    Jth  =  th1(k)*J(:,:,g,1)  +  th2(k)*J(:,:,g,2);
    Gth  =  th1(k)*G(:,:,g,1)  +  th2(k)*G(:,:,g,2);
    Hth  =  th1(k)*H(:,:,g,1)  +  th2(k)*H(:,:,g,2);
    
      
    %% Feedback and Observer Gains and Matrices
    Athh  =  th1h(k)*A(:,:,g,1)  +  th2h(k)*A(:,:,g,2);
    Bthh  =  th1h(k)*B(:,:,g,1)  +  th2h(k)*B(:,:,g,2);
    Cthh  =  th1h(k)*C(:,:,g,1)  +  th2h(k)*C(:,:,g,2);
    
        
    
    Fthh  =  inv(th1h(k)*S(:,:,g,1) + th2h(k)*S(:,:,g,2))*(th1h(k)*F_(:,:,g,1) + th2h(k)*F_(:,:,g,2));
    Lthh  =  inv(th1h(k)*Z(:,:,g,1) + th2h(k)*Z(:,:,g,2))*(th1h(k)*L_(:,:,g,1) + th2h(k)*L_(:,:,g,2));
    
   
    
    u(k)  =  Fthh*xh(:,k);
    %% Observer
    y(:,k)    = Cth*x(:,k)   + Dth*wk;
    xh(:,k+1) = Athh*xh(:,k) + Bthh*u(k) + Lthh*(y(:,k) - Cthh*xh(:,k));
    
    %% System
    x(:,k+1)  = Ath*x(:,k) + Bth*u(k) + Eth*wk;
    z(k)      = Gth*x(:,k) + Hth*u(k) + Jth*wk;
    
    Zsum(k+1) = Zsum(k) + z(k)'*z(k);
    Wsum(k+1) = Wsum(k) + wk'*wk;
    Hinf(k+1) = sqrt(Zsum(k+1)/(gama^2*Wsum(k+1)));

    rho(k+1)  = randsample(1:s,1,true,Pi(g,:));
end

fontsize  = 14;
linewidth = 2;

clf(figure(1)); axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,x(1,:),'linewidth',linewidth); hold on;
stairs(0:N-1,xh(1,:),'linewidth',linewidth);
set(gca,'fontsize',fontsize);
axis([0 N-1 -0.3 0.6]);
legend('$x_{1,k}$','$\hat x_{1,k}$','fontsize',fontsize+6,'interpreter','latex')
grid on;

clf(figure(2));  axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,x(2,:),'linewidth',linewidth); hold on;
stairs(0:N-1,xh(2,:),'linewidth',linewidth);
set(gca,'fontsize',fontsize);
axis([0 N-1 -0.6 0.6]);
legend('$x_{2,k}$','$\hat x_{2,k}$','fontsize',fontsize+6,'interpreter','latex')
grid on;

clf(figure(4)); axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,u,'linewidth',linewidth);
legend('$u(k)$','fontsize',fontsize+6,'interpreter','latex');
set(gca,'fontsize',fontsize);
axis([0 N-1 -0.05 0.01]);
grid on;

clf(figure(5)); axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,Hinf,'linewidth',linewidth);
set(gca,'fontsize',fontsize);
axis([0 N-1 0 1]);
legend('$\frac{\sum_{k=0}^{T}\Vert z_k \Vert_2^2}{\gamma_{min}^2 \sum_{k=0}^{T}\Vert w_k \Vert_2^2}$'...
    ,'fontsize',fontsize+6,'interpreter','latex','location','best');
grid on;

clf(figure(6)); axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,rho,'linewidth',linewidth);
set(gca,'fontsize',fontsize);
legend('$\rho_k$','fontsize',fontsize+6,'interpreter','latex');
axis([0 N-1 0.8 3.2]);
yticks([1 2 3]);
grid on


clf(figure(7)); axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,th1-th1h,'linewidth',linewidth);
set(gca,'fontsize',fontsize);
legend('$|\theta_{i,k} - \hat\theta_{i,k}|$','fontsize',fontsize+6,'interpreter','latex');
grid on


end




