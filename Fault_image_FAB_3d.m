% clear all
% clc
% [Scube,SegyTraceHeaders,SegyHeader]=ReadSegy('Eigenstructure Coherence (3).sgy');
%  Scube=reshape(Scube,[751,440,201]);
%  for it=1:17
%      for icdp=1:17
%         for iline=1:17
%             delt_t=abs(it-9);
%             delt_cdp=abs(icdp-9);
%             delt_line=abs(iline-9);
%             tmp=delt_cdp^2/0.25+delt_t^2/36+delt_line^2/0.25;
%             gau_kern(it,icdp,iline)=exp(-tmp);
%         end
%     end
% end
% figure,imagesc(Scube(:,:,10));colormap(gray)
% Scube_g=convn(Scube,gau_kern,'same');
Scube_g=Scube;
% 
Scube_g=Scube_g/max(Scube_g(:));
figure,imagesc(Scube_g(:,:,10));colormap(gray)
Scube_em=1./(1+exp(-10*(Scube_g-0.8)));
figure,imagesc(Scube_em(:,:,10));colormap(gray)
[nmd1,nmd2,nmd3,dip,dist]=pca_fault_3d(Scube_em,5,5,5);
plane_pro=(nmd2-nmd3)./nmd1;
dip=sin(dip*pi/180);
f_ppro=1./(1+exp(-50*(plane_pro-0.1)));
f_dip=1./(1+exp(-200*(dip-0.98)));
nmd=f_ppro.*f_dip;


kf=0.3;
kb=0.7;
width=0.3;
coff=1./(1+(nmd/kf).^4)-1.5./(1+((nmd-kb)/width).^4);
coffx=derivatives(coff,'x');
coffy=derivatives(coff,'y');
coffz=derivatives(coff,'z');

for i_iter=1:10
    ux=derivatives(Scube_g,'x'); uy=derivatives(Scube_g,'y');
    uz=derivatives(Scube_g,'z');
    uxx=derivatives(ux,'x');uyy=derivatives(uy,'y');
    uzz=derivatives(uz,'z');
    Scube_g=Scube_g+(1/7)*(coff.*(uxx+uyy+uzz)+coffx.*ux+coffy.*uy+coffz.*uz);
    Scube_g=(Scube_g>1).*ones(size(Scube_g))+((Scube_g>0)&(Scube_g<1)).*Scube_g;
end

[output]=Directional_smooth_3d(Scube_g,5,5,5);