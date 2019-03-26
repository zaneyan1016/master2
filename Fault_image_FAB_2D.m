clear all
clc

Coh=ReadSegy('inline20.sgy');



C1=Coh;
C1=C1/max(C1(:));
C1em=1./(1+exp(-10*(C1-0.85)));
[nmd1,nmd2,dip,dist]=pca_fault_2d(C1em,5,5);
line_pro=(nmd1-nmd2)./(nmd1+eps);
f_pro=1./(1+exp(-40*(line_pro-graythresh(line_pro))));
dip=sin(dip*pi/180);
f_dip=1./(1+exp(-20*(dip-0.85)));
f_dist=1./(1+exp(-20*(-dist+1.5)));
nmd=f_dip.*f_pro;

kf=0.3;
kb=0.7;
width=0.3;

coff=1./(1+(nmd/kf).^4)-1.0./(1+((nmd-kb)/width).^4);
coffx=derivatives(coff,'x'); coffy=derivatives(coff,'y');
for i_iter=1:15    
    ux=derivatives(C1,'x'); uy=derivatives(C1,'y');
    uxx=derivatives(ux,'x');uyy=derivatives(uy,'y');
    C1=C1+(1/7)*(coff.*(uxx+uyy)+coffx.*ux+coffy.*uy);
    C1=(C1>1).*ones(size(C1))+((C1>0)&(C1<1)).*C1;
end



C1em=1./(1+exp(-10*(C1-0.8)));
[nmd1,nmd2,dip,dist]=pca_fault_2d(C1em,3,7);
line_pro=(nmd1-nmd2)./(nmd1+eps);
f_pro=1./(1+exp(-40*(line_pro-graythresh(line_pro))));
dip=sin(dip*pi/180);
f_dip=1./(1+exp(-20*(dip-0.85)));
f_dist=1./(1+exp(-20*(-dist+1.5)));
nmd=f_dip.*f_pro;

kf=0.3;
kb=0.7;
width=0.3;
coff=1./(1+(nmd/kf).^4)-1.0./(1+((nmd-kb)/width).^4);
coffx=derivatives(coff,'x'); coffy=derivatives(coff,'y');
for i_iter=1:15
    ux=derivatives(C1,'x'); uy=derivatives(C1,'y');
    uxx=derivatives(ux,'x');uyy=derivatives(uy,'y');
    C1=C1+(1/7)*(coff.*(uxx+uyy)+coffx.*ux+coffy.*uy);
    C1=(C1>1).*ones(size(C1))+((C1>0)&(C1<1)).*C1;
end
figure,imagesc(Coh);colormap(gray);
[output]=Directional_smooth_2d(C1,5,5);
figure,imagesc(output,[0,0.55]);colormap(gray);