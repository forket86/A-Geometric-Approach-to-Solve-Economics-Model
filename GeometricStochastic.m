clear all; clc; close all;

alpha=1/3;
delta=0.1;
F=@(k) k^alpha;
Fk=@(k) alpha*k^(alpha-1);
Fkinv=@(Fk) (Fk/alpha)^(1/(alpha-1));
gamma=2.0;
u=@(c) c^(1-gamma)/(1-gamma);
du=@(c) c^(-gamma);
duinv=@(du) du^(-1/gamma);
beta=0.98;

Kstar=Fkinv(1/beta+delta-1);
K_grd=linspace(Kstar*0.7,Kstar*1.3,10);

Cstar=F(Kstar)+(1-delta)*Kstar-Kstar;
C_grd=linspace(Cstar*0.7,Cstar*1.3,10);
z_grd=[0.95,1.05];
Pz=[0.1 0.9;0.9 0.1];

%% Euler restriction: computer 1
for Knext_iter=1:length(K_grd)
    for KnextnextLOW_iter=1:length(K_grd)
        for KnextnextHIGH_iter=1:length(K_grd)
            for z_iter=1:length(z_grd)
                Knext=K_grd(Knext_iter);
                CnextLOW=z_grd(1)*F(Knext)+(1-delta)*Knext-K_grd(KnextnextLOW_iter);
                CnextHIGH=z_grd(2)*F(Knext)+(1-delta)*Knext-K_grd(KnextnextHIGH_iter);
                C(Knext_iter,KnextnextLOW_iter,KnextnextHIGH_iter,z_iter)=duinv(beta*Pz(z_iter,:)*[du(CnextLOW) du(CnextHIGH)]'*(1+Fk(Knext)-delta));
            end
        end
    end
end
[X,Y,Z,W]=ndgrid(K_grd,K_grd,K_grd,z_grd);
CFE=griddedInterpolant(X,Y,Z,W,C);

%% Budget restriction: computer 2
for K_iter=1:length(K_grd)
    for Knext_iter=1:length(K_grd)
        for z_iter=1:length(z_grd)
            Knext=K_grd(Knext_iter);
            K=K_grd(K_iter);
            CB(K_iter,Knext_iter,z_iter)=z_grd(z_iter)*F(K)+(1-delta)*K-Knext;
        end
    end
end
% figure
% hold on
% surf(K_grd,K_grd,CB(:,:,1),'FaceAlpha',0.4);
% surf(K_grd,K_grd,CB(:,:,2),'FaceAlpha',0.4);
% colormap(bone)
% xlabel('K_{t-1}');
% ylabel('K_{t}');
% zlabel('c_{t}')
% title('Budget')

[X,Y,Z]=ndgrid(K_grd,K_grd,z_grd);
CFB=griddedInterpolant(X,Y,Z,CB);

%% Geometry Intersection: computer 3
opts = optimset('Diagnostics','off', 'Display','off');
KnextFEBintersector=@(K,KnextnextLOW,KnextnextHIGH,z) fzero(@(Knext) CFB(K,Knext,z)-CFE(Knext,KnextnextLOW,KnextnextHIGH,z),Kstar,opts);
for K_iter=1:length(K_grd)
    for KnextnextLOW_iter=1:length(K_grd)
        for KnextnextHIGH_iter=1:length(K_grd)
            for z_iter=1:length(z_grd)
                K=K_grd(K_iter);
                Knextnextlow=K_grd(KnextnextLOW_iter);
                Knextnexthigh=K_grd(KnextnextHIGH_iter);
                z=z_grd(z_iter);
                KnextEB(K_iter,KnextnextLOW_iter,KnextnextHIGH_iter,z_iter)=KnextFEBintersector(K,Knextnextlow,Knextnexthigh,z);
                if isnan(KnextEB(K_iter,KnextnextLOW_iter,KnextnextHIGH_iter,z_iter)) || isreal(KnextEB(K_iter,KnextnextLOW_iter,KnextnextHIGH_iter,z_iter))==false
                    KnextEB(K_iter,KnextnextLOW_iter,KnextnextHIGH_iter,z_iter)=0;
                end
            end
        end
    end
end
[X,Y,Z,W]=ndgrid(K_grd,K_grd,K_grd,z_grd);
KnextFEB=griddedInterpolant(X,Y,Z,W,KnextEB);


%% Graph with Euler-Budget intersection
figure
hold on
surf(z_grd,K_grd,reshape(C(:,1,10,:),length(K_grd),length(z_grd)),'FaceAlpha',0.4);
surf(z_grd,K_grd,reshape(CB(5,:,:),length(K_grd),length(z_grd)),'FaceAlpha',0.4);
z_grd_temp=linspace(z_grd(1),z_grd(end),10);
Knextintersection=reshape(KnextFEB(Kstar.*ones(1,length(z_grd_temp)),K_grd(1).*ones(1,length(z_grd_temp)),K_grd(end).*ones(1,length(z_grd_temp)),z_grd_temp),length(z_grd_temp),1);
Intersect=[z_grd_temp' Knextintersection CFB(Kstar.*ones(length(z_grd_temp),1),Knextintersection,z_grd_temp')];
s=scatter3(Intersect(:,1),Intersect(:,2),Intersect(:,3),'filled')
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'r';
s.MarkerFaceColor = [1 0 0];
hold off
colormap(bone)
shading interp
ylabel('K_{t}');
xlabel('z_t');
zlabel('c_{t}')

%% Restriction to Markov policies: computer 4
phi=zeros(7,1);

Policy=@(K,z,phi) phi(1)*ones(length(z),length(K))+phi(2)*K.*ones(length(z),length(K))+phi(3)*K.^2.*ones(length(z),length(K))+phi(4)*z.*ones(length(z),length(K))+phi(5)*z.^2.*ones(length(z),length(K))+phi(6).*K.*z.*ones(length(z),length(K))+phi(7)*K.^2.*z.^2.*ones(length(z),length(K));
IntersectRestriction=@(phi) Policy(K_grd,z_grd',phi)-KnextFEB(K_grd.*ones(length(z_grd),length(K_grd)),Policy(Policy(K_grd,z_grd',phi),z_grd(1),phi),Policy(Policy(K_grd,z_grd',phi),z_grd(2),phi),ones(length(z_grd),length(K_grd)).*z_grd');
[phi,fval]=fsolve(IntersectRestriction,[Kstar,zeros(1,length(phi)-1)]);

PolicyAnalytical=@(K,z) beta.*alpha.*z.*K.^alpha;
%% Plot & Results
%Check that analytical matches intersection of all restrictions
figure
hold on
fsurf(@(K,z) Policy(K,z,phi),[K_grd(1) K_grd(end) z_grd(1) z_grd(end)]);
fsurf(PolicyAnalytical,[K_grd(1) K_grd(end) z_grd(1) z_grd(end)]); 
%for analytical set delta=1 and log utility
xlabel('K_{t-1}');
ylabel('z_{t}');
zlabel('K_{t}')

%% Graph with the intersection
figure
%subplot(2,1,1)
hold on
surf(K_grd,K_grd,reshape(KnextEB(:,:,1,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
%surf(K_grd,K_grd,reshape(KnextEB(:,:,2,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
surf(K_grd,K_grd,reshape(KnextEB(:,:,3,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
%surf(K_grd,K_grd,reshape(KnextEB(:,:,4,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
surf(K_grd,K_grd,reshape(KnextEB(:,:,5,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
%surf(K_grd,K_grd,reshape(KnextEB(:,:,6,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
surf(K_grd,K_grd,reshape(KnextEB(:,:,7,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
%surf(K_grd,K_grd,reshape(KnextEB(:,:,8,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
%surf(K_grd,K_grd,reshape(KnextEB(:,:,9,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
surf(K_grd,K_grd,reshape(KnextEB(:,:,10,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
K_next_temp=Policy(K_grd,z_grd(1),phi);
K_next_next_temp=Policy(K_next_temp,z_grd(1),phi);
Intersect=[K_grd' K_next_temp' K_next_next_temp'];
s=scatter3(Intersect(:,1),Intersect(:,2),Intersect(:,3),'filled')
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'r';
s.MarkerFaceColor = [1 0 0];
shading interp
colormap(bone)
xlabel('K_{t-1}');
ylabel('K_{t+1}^L');
zlabel('K_{t}')
hold off
% 
%subplot(2,1,2)
figure
hold on
surf(K_grd,K_grd,reshape(KnextEB(:,1,:,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
surf(K_grd,K_grd,reshape(KnextEB(:,3,:,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
surf(K_grd,K_grd,reshape(KnextEB(:,5,:,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
surf(K_grd,K_grd,reshape(KnextEB(:,7,:,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
surf(K_grd,K_grd,reshape(KnextEB(:,10,:,1),length(K_grd),length(K_grd)),'FaceAlpha',0.4);
K_next_temp=Policy(K_grd,z_grd(1),phi);
K_next_next_temp=Policy(K_next_temp,z_grd(2),phi);
Intersect=[K_grd' K_next_next_temp' K_next_temp'];
s=scatter3(Intersect(:,1),Intersect(:,2),Intersect(:,3),'filled')
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'r';
s.MarkerFaceColor = [1 0 0];
shading interp
colormap(bone)
xlabel('K_{t-1}');
ylabel('K_{t+1}^H');
zlabel('K_{t}')
hold off








