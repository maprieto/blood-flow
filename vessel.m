classdef vessel
   %properties
   properties 
      NVAR
      L
      NCELLS
      dx
      xC
      Q
      K
      gamma
      a0
      fnum
      BCL
      BCR
      Qbcs
      Res
   end
   methods
       function v=vessel(NCELLS,L,NVAR,K,gamma,a0,BCL,BCR,Res)
           v.NCELLS=NCELLS;
           v.L=L;
           v.NVAR=NVAR;
           v.K=K;
           v.gamma=gamma;
           v.a0=a0;
           v.dx=L/NCELLS;
           v.xC=zeros(NCELLS,1);
           v.BCL = BCL;
           v.BCR = BCR;
           v.Qbcs = zeros(NVAR,2);
           for i=1:NCELLS
              v.xC(i)=(i-1/2)*v.dx; 
           end
           v.Q=zeros(v.NVAR,v.NCELLS);
           v.Q(1,:)=a0;
           v.Res=Res;
       end
       function plotVessel(v)
           plot(v.xC,v.Q);
       end
   end
end