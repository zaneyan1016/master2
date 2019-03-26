function [nmd1,nmd2,nmd3,dip,dist]=pca_fault_3d(seis_coh,wx,wy,wz)



seis_coh=1-seis_coh;
[nt,ncdp,nline]=size(seis_coh);
wx_oneside=wx;
wy_oneside=wy;
wz_oneside=wz;

new_img=zeros(nt+2*wz_oneside,ncdp+2*wy_oneside,nline+2*wx_oneside);
new_img(wz_oneside+1:nt+wz_oneside,wy_oneside+1:ncdp+wy_oneside,...
    wx_oneside+1:nline+wx_oneside)=seis_coh;
for it=1:wz_oneside
    new_img(it,:,:)=new_img(wz_oneside+1,:,:);
end
for it=nt+wz_oneside+1:nt+2*wz_oneside
    new_img(it,:,:)=new_img(nt+wz_oneside,:,:);
end

for icdp=1:wy_oneside
    new_img(:,icdp,:)=new_img(:,wy_oneside+1,:);
end
for icdp=ncdp+wy_oneside+1:ncdp+2*wy_oneside
    new_img(:,icdp,:)=new_img(:,ncdp+wy_oneside,:);
end

for iline=1:wx_oneside
    new_img(:,:,iline)=new_img(:,:,wx_oneside+1);
end
for iline=nline+wx_oneside+1:nline+2*wx_oneside
    new_img(:,:,iline)=new_img(:,:,nline+wx_oneside);
end

XW=zeros(2*wz_oneside+1,2*wy_oneside+1,2*wx_oneside+1);
ZW=zeros(2*wz_oneside+1,2*wy_oneside+1,2*wx_oneside+1);
YW=zeros(2*wz_oneside+1,2*wy_oneside+1,2*wx_oneside+1);

for icdp=1:2*wy_oneside+1
YW(:,icdp,:)=(icdp-wy_oneside-1);
end
for it=1:2*wz_oneside+1
ZW(it,:,:)=(it-wz_oneside-1);
end
for iline=1:2*wx_oneside+1
XW(:,:,iline)=(iline-wx_oneside-1);
end

parpool(4)
parfor iline=wx_oneside+1:nline+wx_oneside

    iline
    for icdp=wy_oneside+1:ncdp+wy_oneside
        for it=wz_oneside+1:nt+wz_oneside
        a=new_img(it-wz_oneside:it+wz_oneside,...
            icdp-wy_oneside:icdp+wy_oneside,...
            iline-wx_oneside:iline+wx_oneside);
        A=sum(a(:));
        mean_x=sum(sum(sum(a.*XW)))/A;
        mean_z=sum(sum(sum(a.*ZW)))/A;
        mean_y=sum(sum(sum(a.*YW)))/A;
        cxx=sum(sum(sum((XW-mean_x).*(XW-mean_x).*a)))/A;
        cyy=sum(sum(sum((YW-mean_y).*(YW-mean_y).*a)))/A;
        czz=sum(sum(sum((ZW-mean_z).*(ZW-mean_z).*a)))/A;
        cxy=sum(sum(sum((XW-mean_x).*(YW-mean_y).*a)))/A;
        cxz=sum(sum(sum((XW-mean_x).*(ZW-mean_z).*a)))/A;
        cyz=sum(sum(sum((XW-mean_x).*(ZW-mean_z).*a)))/A;

        CM=[cxx,cxy,cxz;cxy,cyy,cyz;cxz,cyz,czz];
        [V,D]=eig(CM);
        distt(it,icdp,iline)=abs(V(1,1)*mean_x+V(2,1)*mean_y+V(3,1)*mean_z);
        D1(it,icdp,iline)=D(3,3);
        D2(it,icdp,iline)=D(2,2);
        D3(it,icdp,iline)=D(1,1);
        dipp(it,icdp,iline)=atan2(sqrt(V(1,1)^2+V(2,1)^2),V(3,1))*180/pi;

        end
    end
end




nmd1=D1(wz_oneside+1:nt+wz_oneside,wy_oneside+1:ncdp+wy_oneside,...
    wx_oneside+1:nline+wx_oneside);
nmd2=D2(wz_oneside+1:nt+wz_oneside,wy_oneside+1:ncdp+wy_oneside,...
    wx_oneside+1:nline+wx_oneside);
nmd3=D3(wz_oneside+1:nt+wz_oneside,wy_oneside+1:ncdp+wy_oneside,...
    wx_oneside+1:nline+wx_oneside);
dip=dipp(wz_oneside+1:nt+wz_oneside,wy_oneside+1:ncdp+wy_oneside,...
    wx_oneside+1:nline+wx_oneside);
dist=distt(wz_oneside+1:nt+wz_oneside,wy_oneside+1:ncdp+wy_oneside,...
    wx_oneside+1:nline+wx_oneside);

