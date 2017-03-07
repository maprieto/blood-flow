classdef model
    properties
        NVAR;
        CFL;
        C;
        nu;
        rho;
        Pv;
    end
    
    methods
        %% PHYSICAL METHOD
        function mod=model(NVAR,CFL)
            mod.NVAR=NVAR;
            mod.CFL=CFL;
            mod.C=vessel.empty(1,0);
            mod.rho=1050.0;
            mod.nu=1.*0.0045/mod.rho;
            mod.Pv=3.2*133.322;
        end
        function mod=add(mod,v)
            mod.C(length(mod.C)+1)=v;
        end
        function F=physicalFlux(mod,Q,i,gamma)
            [a,b]=size(Q);
            F(1,:)=Q(2,:);
            F(2,:)=(Q(2,:).^2./Q(1,:))+gamma*Q(1,:).^(1.5);
            %F(3,:)=Q(3,:).*Q(2,:)./Q(1,:);
        end
        function c=waveSpeed(mod,i,Q)
            c=sqrt(mod.C(i).gamma*3/2*sqrt(Q(1,:)));
        end
        function lambda=eigenvalues(mod,i,Q)
            c=mod.waveSpeed(i,Q);
            
            lambda(1,:)=Q(2,:)./Q(1,:)-c;
            lambda(2,:)=Q(2,:)./Q(1,:)+c;
            %lambda(3,:)=Q(2,:)./Q(1,:);
        end
        function dt=timeStep(mod,i)
            lambda=abs(mod.C(i).Q(2,:)./mod.C(i).Q(1,:)) + mod.waveSpeed(i,mod.C(i).Q);
            maxlambda=max(lambda);
            dt=mod.CFL*mod.C(i).dx/max(maxlambda);
        end
        function a=aFp(mod,K,gamma,a0,p)
           a=(p/K+1.).^2*a0; 
        end
        function p=pFa(mod,K,gamma,a0,Q)
            p=K*((Q(1,:)/a0).^0.5-1.);
        end
        
        %% FLUXES
        function fLF=numFluxLF(mod,i,QL,QR,dt,dx,gamma)
            FL=physicalFlux(mod,QL,i,gamma);
            FR=physicalFlux(mod,QR,i,gamma);
            fLF=0.5*(FL+FR)-0.5*dx/dt*(QR-QL);
        end
        function fLW=numFluxLW(mod,i,QL,QR,dt,dx,gamma)
            FL=physicalFlux(mod,QL,i,gamma);
            FR=physicalFlux(mod,QR,i,gamma);
            QLW=1/2*(QL+QR)-1/2*dt/dx*(FR-FL);
            fLW=physicalFlux(mod,QLW,i,gamma);
        end
        function fLFORCE=numFluxLFORCE(mod,i,QL,QR,dt,dx,gamma)
            fLFORCE=1/2*(numFluxLW(mod,i,QL,QR,dt,dx,gamma)+numFluxLF(mod,i,QL,QR,dt,dx,gamma));
        end
        function fHLLC=numFluxHLLC(mod,i,QL,QR,gamma)      
            aL=QL(1,:);
            aR=QR(1,:);
            uL=QL(2,:)./QL(1,:);
            uR=QR(2,:)./QR(1,:);
%             vL=QL(3,:)./QL(1,:);
%             vR=QR(3,:)./QR(1,:);
            cL=sqrt(mod.C(i).gamma*3/2*sqrt(aL));
            cR=sqrt(mod.C(i).gamma*3/2*sqrt(aR));
            uS=1/2*(uL+uR)-1/2*(cR-cL);
            cS=1/2*(cL+cR)-1/8*(uR-uL);
            sL=min([uL-cL,uS-cS]);
            sR=max([uR+cR,uS+cS]);
            sS=(sL*aR*(uR-sR)-sR*aL*(uL-sL))/(aR*(uR-sR)-aL*(uL-sL));
            aS=aR*(uR-sR)/(sS-sR);
            if (0<=sL)
                fHLLC=physicalFlux(mod,QL,i,gamma);
            elseif (sL<=0 && 0<=sS)
                fHLLC=physicalFlux(mod,QL,i,gamma)+sL*(aS*[1,sS,vL]'-QL);
            elseif (sS<=0 && 0<=sR)
                fHLLC=physicalFlux(mod,QR,i,gamma)+sR*(aS*[1,sS,vR]'-QR);
            elseif (0>=sR)
                fHLLC=physicalFlux(mod,QR,i,gamma);
            end
        end
        function fQvl=numFluxQvl(mod,i,QL,QR,gamma)      
            aL=QL(1,:);
            aR=QR(1,:);
            qL=QL(2,:);
            qR=QR(2,:);
            cL=sqrt(mod.C(i).gamma*3/2*sqrt(aL));
            cR=sqrt(mod.C(i).gamma*3/2*sqrt(aR));
            cM=sqrt(mod.C(i).gamma*3/2*sqrt(0.5*(aL+aR)));
            lambda1M=(qL+qR)/(aL+aR)-cM;
            lambda2M=(qL+qR)/(aL+aR)+cM;
            P=[1 1; lambda1M lambda2M];
            ALambda=[abs(lambda1M) 0;0 abs(lambda2M)];
            MAbsM=P*ALambda*inv(P);
            FL=physicalFlux(mod,QL,i,gamma);
            FR=physicalFlux(mod,QR,i,gamma);
            fQvl=0.5*(FL+FR)-0.5*MAbsM*(QR(1:2,:)-QL(1:2,:));
        end
            function fQvlLLH=numFluxQvlLLH(mod,i,QL,QR,gamma)      
            aL=QL(1,:);
            aR=QR(1,:);
            qL=QL(2,:);
            qR=QR(2,:);
            uL=QL(2,:)./QL(1,:);
            uR=QR(2,:)./QR(1,:);
            cL=sqrt(mod.C(i).gamma*3/2*sqrt(aL));
            cR=sqrt(mod.C(i).gamma*3/2*sqrt(aR));
            cM=sqrt(mod.C(i).gamma*3/2*sqrt(0.5*(aL+aR)));
            lambda1M=(qL+qR)/(aL+aR)-cM;
            lambda2M=(qL+qR)/(aL+aR)+cM;
            lambda1L=qL/aL-cL;
            lambda2L=qL/aL+cR;
            lambda1R=qR/aR-cR;
            lambda2R=qR/aR+cR;
            Regalambda1R=max(abs(lambda1L),abs(lambda1R));
            Regalambda2R=max(abs(lambda2L),abs(lambda2R));
            P=[1 1; lambda1M lambda2M];
            ALambda=[Regalambda1R 0;0 Regalambda2R];
            MAbsM=P*ALambda*inv(P);
            FL=physicalFlux(mod,QL,i,gamma);
            FR=physicalFlux(mod,QR,i,gamma);
            fQvlLLH=0.5*(FL+FR)-0.5*MAbsM*(QR(1:2,:)-QL(1:2,:));
        end
        function fGodEx=numFluxGodunovExact(mod,h,QL,QR,gamma)
            aL=QL(1,:);
            aR=QR(1,:);
            uL=QL(2,:)./QL(1,:);
            uR=QR(2,:)./QR(1,:);
%             vL=QL(3,:)./QL(1,:);
%             vR=QR(3,:)./QR(1,:);
            [aS,uS]=solveERP(mod,h,aL,aR,uL,uR);
            [QGod(1,1),QGod(2,1)]=sampleERP(mod,h,aL,aR,uL,uR,mod.C(h).K,mod.C(h).gamma,aS,uS,0);
            %[QGod(1,1),QGod(2,1),QGod(3,1)]=sampleERP(mod,h,aL,aR,uL,uR,vL,vR,mod.C(h).K,mod.C(h).gamma,aS,uS,0);
            QGod(2,1)=QGod(2,1)*QGod(1,1);
            %QGod(3,1)=QGod(3,1)*QGod(1,1);
            fGodEx=physicalFlux(mod,QGod,h,gamma);
        end
        
        %% EVOLVE FUNCTION
        function Q=evolve(mod,h,dt,dx,gamma)
            S=zeros(2,mod.C(h).NCELLS);
            %S=zeros(3,mod.C(h).NCELLS);
            fNum=zeros(2,mod.C(h).NCELLS+1);
            %fNum=zeros(3,mod.C(h).NCELLS+1);
            for i=1:(mod.C(h).NCELLS+1)
                if i==1
                    QL=mod.C(h).Qbcs(:,1);
                    QR=mod.C(h).Q(:,i);
                elseif i==(mod.C(h).NCELLS+1)
                    QL=mod.C(h).Q(:,i-1);
                    QR=mod.C(h).Qbcs(:,2);
                else
                    QL=mod.C(h).Q(:,i-1);
                    QR=mod.C(h).Q(:,i);
                end
                %fNum(:,i)=numFluxHLLC(mod,h,QL,QR,gamma);
                fNum(:,i)=numFluxGodunovExact(mod,h,QL,QR,gamma); % P
                %fNum(:,i)=numFluxLFORCE(mod,h,QL,QR,dt,dx,gamma);
                %fNum(:,i)=numFluxLW(mod,h,QL,QR,dt,dx,gamma);
                %fNum(:,i)=numFluxLF(mod,h,QL,QR,dt,dx,gamma); % P
                %fNum(:,i)=numFluxQvl(mod,h,QL,QR,gamma); % P
                %fNum(:,i)=numFluxQvlLLH(mod,h,QL,QR,gamma); % P
            end
            Q=mod.C(h).Q-dt/dx*(fNum(:,2:(mod.C(h).NCELLS+1))-fNum(:,1:mod.C(h).NCELLS))+dt*S(:,1:mod.C(h).NCELLS);
        end
        
        %% BOUNDARY
        function [Qbcs,BCL,BCR] = boundaryConditions(mod,h,Q,NCELLS,BCL,BCR,Qbcs,aL,aR,time)
            if BCL==1 % Transparent boundary conditions
                Qbcs(:,1) = Q(:,1);
            elseif BCL==2 % Reflexive boundary conditions
                Qbcs(1,1) = Q(1,1);
                Qbcs(2,1) = -Q(2,1);
%                 Qbcs(3,1) = Q(3,1);
            elseif BCL==3 % Imposed boundary conditions
                Qbcs(1,1) = aL;
                Qbcs(2,1) = UBC(mod,h,aL,-1);
%                 Qbcs(3,1) = Q(3,1);
            elseif BCL==4
                TT = 0.01; % s
                T = 0.827;
                qMax = 100.0*1e-6; % m^3/s
                PI=pi;
                t=time;
                if time <= T*1.2
                    qFixed = qMax*max(sin(pi*time/TT),0.);
                    %qFixed=(3.1199+7.7982*sin(2*PI*t/T+0.5769)+4.1228*sin(4*PI*t/T-0.8738)-1.0611*sin(6*PI*t/T+0.7240)+0.7605*sin(8*PI*t/T-0.6387)-0.9148*sin(10*PI*t/T+1.1598)+0.4924*sin(12*PI*t/T-1.0905)-0.5580*sin(14*PI*t/T+1.042)+0.3280*sin(16*PI*t/T-0.5570)-0.3941*sin(18*PI*t/T+1.2685)+0.2833*sin(20*PI*t/T+0.6702)+0.2272*sin(22*PI*t/T-1.4983)+0.2249*sin(24*PI*t/T+0.9924)+0.2589*sin(26*PI*t/T-1.5616)-0.1460*sin(28*PI*t/T-1.3106)+0.2141*sin(30*PI*t/T-1.1306)-0.1253*sin(32*PI*t/T+0.1552)+0.1321*sin(34*PI*t/T-1.5595)-0.1399*sin(36*PI*t/T+0.4223)-0.0324*sin(38*PI*t/T+0.7811)-0.1211*sin(40*PI*t/T+1.0729))/1000/60;
                else
                    %qFixed = 0.;
                    %mod.C(h).BCL = 2;
                end
                qFixed=(3.1199+7.7982*sin(2*PI*t/T+0.5769)+4.1228*sin(4*PI*t/T-0.8738)-1.0611*sin(6*PI*t/T+0.7240)+0.7605*sin(8*PI*t/T-0.6387)-0.9148*sin(10*PI*t/T+1.1598)+0.4924*sin(12*PI*t/T-1.0905)-0.5580*sin(14*PI*t/T+1.042)+0.3280*sin(16*PI*t/T-0.5570)-0.3941*sin(18*PI*t/T+1.2685)+0.2833*sin(20*PI*t/T+0.6702)+0.2272*sin(22*PI*t/T-1.4983)+0.2249*sin(24*PI*t/T+0.9924)+0.2589*sin(26*PI*t/T-1.5616)-0.1460*sin(28*PI*t/T-1.3106)+0.2141*sin(30*PI*t/T-1.1306)-0.1253*sin(32*PI*t/T+0.1552)+0.1321*sin(34*PI*t/T-1.5595)-0.1399*sin(36*PI*t/T+0.4223)-0.0324*sin(38*PI*t/T+0.7811)-0.1211*sin(40*PI*t/T+1.0729))/1000/60;
                Qbcs(2,1) = qFixed;
                Qbcs(1,1) = AQBC(mod,h,qFixed,-1);
%                 Qbcs(3,1) = Q(3,1);
            elseif BCL==5
                TT=0.010;
                pMax=30.0*133.32;
                if time<=TT*1.2
                   pFixed=80.*133.32+pMax*max(sin(pi*time/TT),0.);
                   aFixed=aFp(mod,mod.C(h).K,mod.C(h).gamma,mod.C(h).a0,pFixed);
                   Qbcs(1,1)=aFixed;
                   %rar=RiemannInv(mod,K,gamma,Q(:,1),-1);
                   Qbcs(2,1)=UBC(mod,h,aFixed,-1);%inflowFixedP(mod,K,gamma,aFixed,Q(:,1),rar);
 %                  Qbcs(3,1)=Q(3,1);
                else
                   Qbcs(1,1) = Q(1,1);
                   Qbcs(2,1) = 0;
 %                  Qbcs(3,1) = Q(3,1); 
                   mod.C(h).BCL=2;
                end
            end
            if BCR==1 % Transparent boundary conditions
                Qbcs(:,2) = Q(:,NCELLS);
            elseif BCR==2 % Reflexive boundary conditions
                Qbcs(1,2) = Q(1,NCELLS);
                Qbcs(2,2) = -Q(2,NCELLS);
 %               Qbcs(3,2) = Q(3,NCELLS);
            elseif BCR==3 % Reflexive boundary conditions
                Qbcs(1,2) = aR;
                Qbcs(2,2) = UBC(mod,h,aR,1);
 %               Qbcs(3,2) = Q(3,NCELLS);
            elseif BCR==4
                TT = 0.01; % s
                qMax = 100.0*1e-6; % m^3/s
                if time <= TT*1.2
                    qFixed = qMax*max(sin(pi*time/TT),0.);
                else
                    qFixed = 0.;
                    mod.C(h).BCR = 2;
                end
                
                Qbcs(2,2) = qFixed;
                Qbcs(1,2) = ABC(mod,h,qFixed,1);
   %             Qbcs(3,2) = Q(3,2);
            elseif BCR==5
                TT=0.01;
                pMax=10.0*133.32;
                if time<=TT*1.2
                   pFixed=80.*133.32+pMax*max(sin(pi*time/TT),0.);
                   aFixed=aFp(mod,mod.C(h).K,mod.C(h).gamma,mod.C(h).a0,pFixed);
                   Qbcs(1,2)=aFixed;
                   %rar=RiemannInv(mod,K,gamma,Q(:,1),1);
                   Qbcs(2,2)=UBC(mod,h,aFixed,1);
   %                Qbcs(3,2)=Q(3,1);
                else
                   Qbcs(1,2) = Q(1,1);
                   Qbcs(2,2) = 0;
   %                Qbcs(3,2) = Q(3,1); 
                   mod.C(h).BCR=2;
                end
            elseif BCR==6
                n=mod.C(h).NCELLS;
                cost=mod.RiemannInv(h,1);
                options = optimoptions('fsolve','Display','off');
                x0=[mod.C(h).Q(1,n)];
                p=@(a)pFa(mod,mod.C(h).K,mod.C(h).gamma,mod.C(h).a0,a);
                c=@(a)waveSpeed(mod,h,a);
                f=@(a)cost*a*mod.C(h).Res-4*c(a)*a*mod.C(h).Res-p(a)+mod.Pv;
                x=fsolve(@(a) f(a),x0,options);
                Qbcs(1,2) = x;
                Qbcs(2,2) = (p(x)-mod.Pv)/(mod.C(h).Res); %AQBC(mod,h,x,1);
                %Qbcs(2,2) = UBC(mod,h,x,1);
  %              Qbcs(3,2) = 0;
            end
        end

        %% EXACT SOLVER
        function [arar,urar]=rarL(mod,i,uL,cL,xt)
            urar=1/5*(uL+4*cL+4*xt);
            crar=1/5*(uL+4*cL-xt);
            arar=crar^4*4/9/mod.C(i).gamma^2;
            %arar=(2/3/mod.C(i).gamma*(1.5.*(4*cR*uR-xt))^2)^2;
        end
        function [arar,urar]=rarR(mod,i,uR,cR,xt)
            urar=1/5*(uR-4*cR+4*xt);
            %disp(urar);
            crar=1/5*(-uR+4*cR+xt);
            
            arar=crar^4*4/9/(mod.C(i).gamma^2);
            %arar=(2/3/mod.C(i).gamma*(1.5.*(4*cR*uR-xt))^2)^2;
        end
        function f=fK(mod,i,gamma,aK,aS,cK)
            %cS=waveSpeed(mod,i,aS);
            cS=sqrt(gamma*3/2*sqrt(aS));
            if (aS<=aK)
                f=4*(cS-cK);
            else
                f=sqrt(gamma*(aS-aK)*(aS^(1.5)-aK^(1.5))/(aK*aS));
            end
        end
        function [a,u]=solveERP(mod,i,aL,aR,uL,uR)
            cL=waveSpeed(mod,i,aL);
            cR=waveSpeed(mod,i,aR);
            tot=1e-6;
            eps=1e-12;
            A=1/2*(aL+aR);
            err=1;
            gamma=mod.C(i).gamma;
            while (err>tot)
                fp=(fK(mod,i,gamma,aL,A+eps,cL)+fK(mod,i,gamma,aR,A+eps,cR));
                fm=(fK(mod,i,gamma,aL,A-eps,cL)+fK(mod,i,gamma,aR,A-eps,cR));
                fder=(fp-fm)/(2*eps);
                Aold=A;%AK
                A=A-(fK(mod,i,gamma,aL,A,cL)+fK(mod,i,gamma,aR,A,cR)+uR-uL)/fder;%Ak+1
                err=abs(A-Aold)/((A+Aold)/2);
            end
            a=A;
            u=1/2*(uL+uR)+1/2*(fK(mod,i,gamma,aR,A,cR)-fK(mod,i,gamma,aL,A,cL));
            
        end
        function m = mK(mod,gamma,aK,aS)
            m=sqrt(gamma*aS*aK*(aS^(1.5)-aK^(1.5))/(aS-aK));
        end
        function [a,u] = sampleERP(mod,i,aL,aR,uL,uR,K,gamma,aS,uS,xt)
       %function [a,u,v] = sampleERP(mod,i,aL,aR,uL,uR,vL,vR,K,gamma,aS,uS,xt)
            %cS=waveSpeed(mod,i,aS);
            %cL=waveSpeed(mod,i,aL);
            %cR=waveSpeed(mod,i,aR);
            cS=sqrt(gamma*3/2*sqrt(aS));
            cL=sqrt(gamma*3/2*sqrt(aL));
            cR=sqrt(gamma*3/2*sqrt(aR));
            mL=sqrt(gamma*aS*aL*(aS^(1.5)-aL^(1.5))/(aS-aL));
            mR=sqrt(gamma*aS*aR*(aS^(1.5)-aR^(1.5))/(aS-aR));
            sL=uL-mL/aL;
            sR=uR+mR/aR;
            if (xt<uS)
                if (aS<=aL)
                    if (xt<(uL-cL))
                        a=aL;
                        u=uL;
                    elseif (xt>=(uL-cL) && xt<(uS-cS))
                        [a,u]=rarL(mod,i,uL,cL,xt);
                        
                    elseif (xt>=(uS-cS))
                        a=aS;
                        u=uS;
                        
                    end
                else
                    if (xt<sL)
                        a=aL;
                        u=uL;
                    else
                        a=aS;
                        u=uS;
                        
                    end
                end
            else
                
                if (aS<=aR)
                    if (xt<(uS+cS))
                        a=aS;
                        u=uS;
                    elseif (xt>=(uS+cS) && xt<(uR+cR))
                        [a,u]=rarR(mod,i,uR,cR,xt);
                    elseif (xt>=(uR+cR))
                        a=aR;
                        u=uR;
                    end
                else
                    if (xt<sR)
                        a=aS;
                        u=uS;
                    else
                        a=aR;
                        u=uR;
                    end
                end
            end
            
            
        end
        function Q = exactSampleERP(mod,aL,aR,uL,uR,h,aS,uS,time,gate)
       %function Q = exactSampleERP(mod,aL,aR,uL,uR,vL,vR,h,aS,uS,time,gate)
            NCELLS=mod.C(h).NCELLS;
            xC=mod.C(h).xC;
            K=mod.C(h).K;
            gamma=mod.C(h).gamma;
            Q=zeros(2,NCELLS);
            for i = 1:NCELLS
                [Q(1,i),Q(2,i)] = mod.sampleERP(h,aL,aR,uL,uR,K,gamma,aS,uS,(xC(i)-gate)/(time+1e-10));
                %[Q(1,i),Q(2,i),Q(3,i)] = mod.sampleERP(h,aL,aR,uL,uR,vL,vR,K,gamma,aS,uS,(xC(i)-gate)/(time+1e-10));
            end
        end
    end
end
