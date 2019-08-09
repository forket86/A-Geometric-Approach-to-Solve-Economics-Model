clear all; clc;

alpha=1/2;
delta=0.07;
F=@(k) k^alpha;
Fk=@(k) alpha*k^(alpha-1);
Fkinv=@(Fk) (Fk/alpha)^(1/(alpha-1));
gamma=2;
u=@(c) c^(1-gamma)/(1-gamma);
du=@(c) c^(-gamma);
duinv=@(du) du^(-1/gamma);
beta=0.98;

Kstar=Fkinv(1/beta+delta-1);
K_grd=linspace(Kstar*0.7,Kstar*1.3,10);

for Knext_iter=1:length(K_grd)
    for K_iter=1:length(K_grd)
        K=K_grd(K_iter);
        Knext=K_grd(Knext_iter);
        C=F(K)+(1-delta)*K-Knext;
        
        Knextnext(K_iter,Knext_iter)=-(duinv(du(C)/(beta*(1+Fk(Knext)-delta)))-(F(Knext)+(1-delta)*Knext));
        %% Check
        %Cnext=F(Knext)+(1-delta)*Knext-Knextnext(K_iter,Knext_iter);
        %check(K_iter,Knext_iter)=du(C)-beta*(du(Cnext)*(1+Fk(Knext)-delta));
    end
end

[X,Y]=ndgrid(K_grd,K_grd);
KnextnextFE=griddedInterpolant(X,Y,Knextnext);

hold on
surf(K_grd,K_grd,Knextnext,'FaceAlpha',0.4);
shading interp
colormap(bone)
%title('Geometric Restriction and Policy')
xlabel('K_{t-1}');
ylabel('K_{t}');
zlabel('K_{t+1}');


phi=zeros(3,1);
Policy=@(K,phi) phi(1)*ones(1,length(K))+phi(2)*K+phi(3)*K.^2;

IntersectRestriction=@(phi) Policy(Policy(K_grd,phi),phi)-KnextnextFE(K_grd,Policy(K_grd,phi));

[phi,fval]=fsolve(IntersectRestriction,[Kstar,zeros(1,length(phi)-1)]);


KnextI=Policy(K_grd,phi);
KnextnextI=Policy(KnextI,phi);

Intersect=[K_grd' KnextI' KnextnextI'];
s=scatter3(Intersect(:,1),Intersect(:,2),Intersect(:,3),'filled')
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'r';
s.MarkerFaceColor = [1 0 0];
hold off







