function blood_flow_compute(aL_aux,aR_aux,uL,uR,NCELLS,CFL,handles)

%% DEFINE PARAMETERS FOR VESSELSdefine parameters for vessels
L=0.2; %m
E=4.0e5; %pascal
h=1.1e-3; %m
r0=9.99e-3; %m
a0=pi*r0^2; %m^2
rho=1050; %kg/m^3
gate=0.1; %m
%NCELLS=80;
% Arterias
m=0.5;
K=sqrt(pi)/(1-0.5^2)*E*h/sqrt(a0);
gamma=m*K/rho/(m+1)/a0^m;

%% DEFINE THE VESSELS
v=vessel(NCELLS,L,2,K,gamma,a0,1,1,10);
vessEx = vessel(NCELLS,L,2,K,gamma,a0,1,1,10);
% v=vessel(NCELLS,L,3,K,gamma,a0,1,1,10);
% vessEx = vessel(500,L,3,K,gamma,a0,1,1,10);

%% DEFINE THE MODEL
%CFL=0.7;
mod=model(2,CFL);
mod=mod.add(v);
mod=mod.add(vessEx);

%% INITIAL CONDITION
tEnd = 0.014;
gate = 0.1;
aL = aL_aux*v.a0;
aR = aR_aux*v.a0;

% RP = 1 ;
% if RP == 1 % Rar-shock
%     aL = 2.2*v.a0;
%     aR = 1.1*v.a0;
% %     vL=1;
% %     vR=0.0;
%     uL = 0.;
%     uR = 0.;
%     tEnd = 0.014;
%     gate = 0.1;
% elseif RP == 2 % Rar-rar
%     aL = 1.6*v.a0;
%     aR = 1.6*v.a0;
% %     vL=1;
% %     vR=0.;
%     uL = -1.5;
%     uR = 1.5;
%     tEnd = 0.01;
%     gate = 0.1;
% elseif RP == 3 % Shock-shock
%     aL = 1.6*v.a0;
%     aR = 1.6*v.a0;
% %     vL=1.;
% %     vR=0.;
%     uL = 1.5;
%     uR = -1.5;
%     tEnd = 0.014;
%     gate = 0.1;
% elseif RP == 4 % Shock
% %   AL      = 3.E-4;   %area left
%     aL = 1.6*v.a0;
%     aR = 1.6*v.a0;
% %     vL=1.;
% %     vR=0.;
%     uL = -2.6575E-5;
%     uR = 6.1230E-6;
%     tEnd = 0.014;
%     gate = 0.1;
% elseif RP == 5 % Sonic point
% %   AL      = 3.E-4;   %area left
%     aL = 10.*v.a0;
%     aR = 1.*v.a0;
% %     vL=1.;
% %     vR=0.;
%     uL = 0;
%     uR = 0;
%     tEnd = 0.007;
%     gate = 0.1;
% end

%timing
time=0.;
dt=0.;
nMax=100000; %maximum number of time iterations

%% EXACT VESSEL
[aS,uS] = mod.solveERP(1,aL,aR,uL,uR);
exact = zeros(mod.C(2).NCELLS,3);
exact(:,1) = mod.C(2).xC';
mod.C(2).Q = mod.exactSampleERP(aL,aR,uL,uR,2,aS,uS,time,gate);
%mod.C(2).Q = mod.exactSampleERP(aL,aR,uL,uR,vL,vR,2,aS,uS,time,gate);
exact(:,2) = (mod.C(2).Q(1,:)/mod.C(2).a0)';
exact(:,3) = mod.C(2).Q(2,:)';
% exact(:,4) = mod.C(2).Q(3,:)';

%% INITIAL CONDITION
mod.C(1).Q(1,:)=(mod.C(1).xC<=gate)*aL+(mod.C(1).xC>gate)*aR;
mod.C(1).Q(2,:)=(mod.C(1).xC<=gate)*aL*uL+(mod.C(1).xC>gate)*aR*uR;
% mod.C(1).Q(3,:)=(mod.C(1).xC<=gate)*aL*vL+(mod.C(1).xC>gate)*aR*vR;

%% FIGURE
xplot = mod.C(1).xC/mod.C(1).L;
aplot = mod.C(1).Q(1,:)'/mod.C(1).a0;
uplot = mod.C(1).Q(2,:)'./mod.C(1).Q(1,:)';
ruplot = sqrt(mod.C(1).Q(1,:)'./mod.C(1).a0./pi)';
rlplot = -sqrt(mod.C(1).Q(1,:)'./mod.C(1).a0./pi)';
% vplot = mod.C(1).Q(3,:)'./mod.C(1).Q(1,:)';
eplot = mod.eigenvalues(1,mod.C(1).Q);

%% EIGENVALUE
cS=sqrt(mod.C(1).gamma*3/2*sqrt(aS));
cL=sqrt(mod.C(1).gamma*3/2*sqrt(aL));
cR=sqrt(mod.C(1).gamma*3/2*sqrt(aR));
T=linspace(0,tEnd,NCELLS);
sL=uL-mK(mod,gamma,aL,aS)/aL;
sR=uR+mK(mod,gamma,aR,aS)/aR;
uLplot =T.*(sL);
uRplot =T.*(sR);
uCplot =T.*(uS);
Uplot(1,:)=uLplot;
Uplot(2,:)=uRplot;
%Uplot(3,:)=uCplot;
if (aS<=aR)
    [a,b]=size(Uplot);
    uRplot =T.*(uR+cR);
    uR2plot=T.*(uS+cS);
    Uplot(a+1,:)=uR2plot;
end
if (aS<=aL)
    [a,b]=size(Uplot);
    uLplot =T.*(uS-cS);
    uL2plot=T.*(uL-cL);
    Uplot(a+1,:)=uL2plot;
end
Uplot(1,:)=uLplot;
Uplot(2,:)=uRplot;
%Uplot(3,:)=uCplot;
timeplot=zeros(1,NCELLS);
timeplot(:)=time;

%% PLOT
if(get(handles.checkbox1,'Value')==0)
    %%%%%%%%%%%%  2D plot %%%%%%%%%%%%%%%%
    cla(handles.axes9)
    set(handles.axes9,'Visible','off')
    set(handles.axes1,'Visible','on')
    set(handles.axes2,'Visible','on')
    set(handles.axes3,'Visible','on')
    set(handles.axes4,'Visible','on')
    set(handles.text9,'Visible','on')
    
    cla(handles.axes1)
    axes(handles.axes1)
    hold(handles.axes1,'on')
    aLine = plot(xplot, aplot,'o', 'YDataSource', 'aplot', 'XDataSource', 'xplot');
    aExactLine = plot(exact(:,1)/mod.C(2).L, exact(:,2),'r-', 'XDataSource', ...
        'exact(:,1)/mod.C(2).L', 'YDataSource', 'exact(:,2)');
    hold(handles.axes1,'off')
    %xlabel(handles.axes1,'x/L')
    ylabel(handles.axes1,'A/A0')
    
    % Grafic Representation
    cla(handles.axes2)
    axes(handles.axes2)
    hold(handles.axes2,'on')
    ruLine = plot(xplot, ruplot,'o', 'YDataSource', 'ruplot', 'XDataSource', 'xplot');
    rlLine = plot(xplot, rlplot,'o', 'YDataSource', 'rlplot', 'XDataSource', 'xplot');
    ruExactLine = plot(exact(:,1)/mod.C(2).L, sqrt(exact(:,2)/pi),'r-', 'XDataSource', ...
        'exact(:,1)/mod.C(2).L', 'YDataSource', 'sqrt(exact(:,2)/pi)');
    rlExactLine = plot(exact(:,1)/mod.C(2).L, -sqrt(exact(:,2)/pi),'r-', 'XDataSource', ...
        'exact(:,1)/mod.C(2).L', 'YDataSource', '-sqrt(exact(:,2)/pi)');
    hold(handles.axes2,'off')
    %xlabel('x [m]')
    ylabel(handles.axes2,'D - vessel Diameter [m]')
    
    cla(handles.axes3)
    axes(handles.axes3)
    hold(handles.axes3,'on')
    uLine = plot(xplot, uplot,'o', 'YDataSource', 'uplot', 'XDataSource', 'xplot');
    uExactLine = plot(exact(:,1)/mod.C(2).L, exact(:,3),'r-', 'XDataSource', 'exact(:,1)/mod.C(2).L', 'YDataSource', 'exact(:,3)');
    hold(handles.axes3,'off')
    %xlabel('x/L')
    ylabel(handles.axes3,'u [m/s]')
    
    % subplot(3,1,3)
    % %subplot(4,1,3)
    % vLine = plot(xplot, vplot,'o', 'YDataSource', 'vplot', 'XDataSource', 'xplot');
    % hold on;
    % vExactLine = plot(exact(:,1)/mod.C(1).L, exact(:,4),'r-', 'XDataSource', ...
    %     'exact(:,1)/mod.C(1).L', 'YDataSource', 'exact(:,4)');
    % hold on;
    % axis([0,1,min(vL,vR)-0.3,max(vL,vR)+0.3])
    % xlabel('x/L')
    % ylabel('v')
    
    cla(handles.axes4)
    axes(handles.axes4)
    hold(handles.axes4,'on')
    eLine= plot((Uplot+gate)./mod.C(2).L,T);
    timeLine=plot(linspace(0,1,NCELLS),timeplot,'YdataSource','timeplot','XdataSource','linspace(0,1,NCELLS)');
    axis(handles.axes4,[0,1,0,tEnd])
    hold(handles.axes4,'off')
    %legend('A/A0','U','C');
    
else
    %%%%%%% 3D plot %%%%%%%%%%%%%%
    cla(handles.axes1)
    cla(handles.axes2)
    cla(handles.axes3)
    cla(handles.axes4)
    set(handles.axes1,'Visible','off')
    set(handles.axes2,'Visible','off')
    set(handles.axes3,'Visible','off')
    set(handles.axes4,'Visible','off')
    set(handles.text9,'Visible','off')
    set(handles.text8,'string','Time: 0 s');
    set(handles.axes9,'Visible','on')
    cla(handles.axes9)
    axes(handles.axes9)
    % Radius profile and surface plot
    r=ruplot; % Radius profile in 1D
    n_theta=100;
    len=1;
    [X,Y,Z]=cylinder(handles.axes9,r,n_theta);
    hsurf=surf(handles.axes9,X,len*Z,Y);
    
    % Setup visualization
    lmax=max(abs(sqrt([aL_aux,aR_aux]/pi)));
    axis(handles.axes9,[-lmax,lmax,0,1,-lmax,lmax])
    shading(handles.axes9,'interp')
    %xlabel('x'),xlabel('y'),xlabel('z')
    colormap(handles.axes9,'copper'), shading(handles.axes9,'interp')
    light('Position',[20,2,-20]),lighting phong, material([0.4,0.6,0.5,30])
    set(hsurf,'FaceColor','y','BackFaceLighting','lit','FaceAlpha',0.6)
    view(handles.axes9,[-60,20]),set(handles.axes9,'CameraViewAngleMode','Manual')
    
    % Velocity field
    v=uplot; % velocity array in 1D
    
    % Add velocity field in inner cylindrical sublayers
    amp=[0.74, 0.55, 0.35, 0.];
    n_theta_plot=12;r_gap=floor(numel(r)/10);
    V=repmat([v(1),v(2:r_gap:end-1)',v(end)],1,n_theta_plot+1);
    X1=cos(repmat(linspace(0,2*pi,n_theta_plot+1)',1,numel([r(1),r(2:r_gap:end-1),r(end)])));
    Y1=sin(repmat(linspace(0,2*pi,n_theta_plot+1)',1,numel([r(1),r(2:r_gap:end-1),r(end)])));
    R=repmat([r(1),r(2:r_gap:end-1),r(end)]',1,n_theta_plot+1);
    Z=repmat(linspace(0,1,numel([r(1),r(2:r_gap:end-1),r(end)]))',1,n_theta_plot+1);
    hold(handles.axes9,'on')
    for j=1:numel(amp)
        %[X,Y,Z]=cylinder(r(1:r_gap:end)*amp(j),n_theta_plot);
        %[x, y, z, u, v, w, cx, cy, cz, s, color, quiv, method, nointerp]
        hcone=coneplot(handles.axes9,amp(j)*R.*X1',len*Z,amp(j)*R.*Y1',0*X1,V,0*Y1, 0.1,'r','nointerp');
        %colormap jet
        %set(hcone,'FaceColor','interp','Cdata',V)
        %set(hcone,'FaceColor','red')
    end
    hold(handles.axes9,'off')
    set(handles.text8,'string',['Time: ',num2str(time,'%5.4f'),' s']);
    drawnow
end

for k=1:nMax
    if time>=tEnd
        break
    end
    
    dt=mod.timeStep(1);             %compute the time step
    dx=mod.C(1).dx;
    
    if time+dt>tEnd                 %ensure exact exit time matching
        dt=tEnd-time;
    end
    %evolve the solution
    [mod.C(1).Qbcs,mod.C(1).BCL,mod.C(1).BCR] = mod.boundaryConditions(1,mod.C(1).Q,mod.C(1).NCELLS,mod.C(1).BCL,mod.C(1).BCR,mod.C(1).Qbcs,1/2*(aL+aR),1/2*(aL+aR));
    mod.C(1).Q = mod.evolve(1,dt,dx,mod.C(1).gamma);
    %update the time
    time=time+dt;
    timeplot(:)=time;
    
    mod.C(2).Q = mod.exactSampleERP(aL,aR,uL,uR,2,aS,uS,time,gate);
    %mod.C(2).Q = mod.exactSampleERP(aL,aR,uL,uR,vL,vR,2,aS,uS,time,gate);%Update Exact solution
    exact(:,2) = (mod.C(2).Q(1,:)/mod.C(2).a0)';
    exact(:,3) = mod.C(2).Q(2,:)';
    %exact(:,4) = mod.C(2).Q(3,:)';
    % Refresh plot
    aplot = mod.C(1).Q(1,:)'/mod.C(1).a0;
    uplot = mod.C(1).Q(2,:)'./mod.C(1).Q(1,:)';
    %rplot = sqrt(mod.C(1).Q(1,:)'./mod.C(1).a0/pi)';
    rlplot = sqrt(aplot/pi)';
    ruplot = -sqrt(aplot/pi)';
    %vplot = mod.C(1).Q(3,:)'./mod.C(1).Q(1,:)';
    
    if(get(handles.checkbox1,'Value')==0)
        %%%%%%%%%%%%  2D plot %%%%%%%%%%%%%%%%
        axes(handles.axes1)
        aLine = plot(xplot, aplot,'o', 'YDataSource', 'aplot', 'XDataSource', 'xplot');
        hold(handles.axes1,'on')
        aExactLine = plot(exact(:,1)/mod.C(2).L, exact(:,2),'r-', 'XDataSource', ...
            'exact(:,1)/mod.C(2).L', 'YDataSource', 'exact(:,2)');
        hold(handles.axes1,'off')
        ylabel(handles.axes1,'A/A0')
        
        axes(handles.axes2)
        ruLine = plot(xplot, ruplot,'o', 'YDataSource', 'ruplot', 'XDataSource', 'xplot');
        hold(handles.axes2,'on')
        rlLine = plot(xplot, rlplot,'o', 'YDataSource', 'rlplot', 'XDataSource', 'xplot');
        ruExactLine = plot(exact(:,1)/mod.C(2).L, sqrt(exact(:,2)/pi),'r-', 'XDataSource', ...
            'exact(:,1)/mod.C(2).L', 'YDataSource', 'sqrt(exact(:,2)/pi)');
        rlExactLine = plot(exact(:,1)/mod.C(2).L, -sqrt(exact(:,2)/pi),'r-', 'XDataSource', ...
            'exact(:,1)/mod.C(2).L', 'YDataSource', '-sqrt(exact(:,2)/pi)');
        hold(handles.axes2,'off')
        ylabel(handles.axes2,'D - vessel Diameter [m]')
        
        axes(handles.axes3)
        uLine = plot(xplot, uplot,'o', 'YDataSource', 'uplot', 'XDataSource', 'xplot');
        hold(handles.axes3,'on')
        uExactLine = plot(exact(:,1)/mod.C(2).L, exact(:,3),'r-', 'XDataSource', 'exact(:,1)/mod.C(2).L', 'YDataSource', 'exact(:,3)');
        hold(handles.axes3,'off')
        ylabel(handles.axes3,'u [m/s]')
        
        %refreshdata(timeLine);
        axes(handles.axes4)
        eLine= plot((Uplot+gate)./mod.C(2).L,T);
        hold(handles.axes4,'on')
        timeLine=plot(linspace(0,1,NCELLS),timeplot,'YdataSource','timeplot','XdataSource','linspace(0,1,NCELLS)');
        axis(handles.axes4,[0,1,0,tEnd])
        hold(handles.axes4,'off')
        
        set(handles.text8,'string',['Time: ',num2str(time,'%5.4f'),' s']);
        drawnow
        
    else
        %%%%%%%%%%%%  3D plot %%%%%%%%%%%%%%%%
        cla(handles.axes9)
        % Radius profile and surface plot
        r=ruplot; % Radius profile in 1D
        n_theta=100;
        len=1;
        [X,Y,Z]=cylinder(handles.axes9,r,n_theta);
        hsurf=surf(handles.axes9,X,len*Z,Y);
        
        % Setup visualization
        lmax=max(abs(sqrt([aL_aux,aR_aux]/pi)));
        axis(handles.axes9,[-lmax,lmax,0,1,-lmax,lmax])
        shading(handles.axes9,'interp')
        %xlabel('x'),xlabel('y'),xlabel('z')
        colormap(handles.axes9,'copper'), shading(handles.axes9,'interp')
        light('Position',[20,2,-20]),lighting phong, material([0.4,0.6,0.5,30])
        set(hsurf,'FaceColor','y','BackFaceLighting','lit','FaceAlpha',0.6)
        view(handles.axes9,[-60,20]),set(handles.axes9,'CameraViewAngleMode','Manual')
        
        % Velocity field
        v=uplot; % velocity array in 1D
        
        % Add velocity field in inner cylindrical sublayers
        amp=[0.74, 0.55, 0.35, 0.];
        n_theta_plot=12;r_gap=floor(numel(r)/10);
        V=repmat([v(1),v(2:r_gap:end-1)',v(end)],1,n_theta_plot+1);
        X1=cos(repmat(linspace(0,2*pi,n_theta_plot+1)',1,numel([r(1),r(2:r_gap:end-1),r(end)])));
        Y1=sin(repmat(linspace(0,2*pi,n_theta_plot+1)',1,numel([r(1),r(2:r_gap:end-1),r(end)])));
        R=repmat([r(1),r(2:r_gap:end-1),r(end)]',1,n_theta_plot+1);
        Z=repmat(linspace(0,1,numel([r(1),r(2:r_gap:end-1),r(end)]))',1,n_theta_plot+1);
        hold(handles.axes9,'on')
        for j=1:numel(amp)
            %[X,Y,Z]=cylinder(r(1:r_gap:end)*amp(j),n_theta_plot);
            %[x, y, z, u, v, w, cx, cy, cz, s, color, quiv, method, nointerp]
            hcone=coneplot(handles.axes9,amp(j)*R.*X1',len*Z,amp(j)*R.*Y1',0*X1,V,0*Y1, 0.1,'r','nointerp');
            %colormap jet
            %set(hcone,'FaceColor','interp','Cdata',V)
            %set(hcone,'FaceColor','red')
        end
        hold(handles.axes9,'off')
        set(handles.text8,'string',['Time: ',num2str(time,'%5.4f'),' s']);
        drawnow
        
    end
    
end
error=zeros(2,NCELLS);
error(1,:)=abs(mod.C(2).Q(1,:)-mod.C(1).Q(1,:));
error(2,:)=abs(mod.C(2).Q(1,:).*mod.C(2).Q(2,:)-mod.C(1).Q(2,:));
l1error=zeros(2,1);
l1error(1)=mod.C(1).dx*sum(error(1,:));
l1error(2)=mod.C(1).dx*sum(error(2,:));
disp(['l1error q_1=A =',num2str(l1error(1))])
disp(['l1error q_2=Au=',num2str(l1error(2))])
%
linferror=zeros(2,1);
linferror(1)=max(error(1,:));
linferror(2)=max(error(2,:));
disp(['linferror(1)=',num2str(linferror(1))])
disp(['linferror(2)=',num2str(linferror(2))])
