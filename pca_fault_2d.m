function [nmd1,nmd2,dip,dist]=pca_fault_2d(seis_coh,wx,wy)

seis_coh=1-seis_coh;

[nz,nx]=size(seis_coh);

wxl_single=wx;
wxl_double=2*wxl_single+1;
wzl_single=wy;
wzl_double=2*wzl_single+1;
new_img=zeros(nz+2*wzl_single,nx+2*wxl_single);
new_img(wzl_single+1:nz+wzl_single,wxl_single+1:nx+wxl_single)=seis_coh;
for iz=1:wzl_single
    new_img(iz,:)=new_img(wzl_single+1,:);
end
for iz=nz+wzl_single+1:nz+2*wzl_single
    new_img(iz,:)=new_img(nz+wzl_single,:);
end

for ix=1:wxl_single
    new_img(:,ix)=new_img(:,wxl_single+1);
end
for ix=nx+wxl_single+1:nx+2*wxl_single
    new_img(:,ix)=new_img(:,nx+wxl_single);
end

XW=zeros(wzl_double,wxl_double);
ZW=zeros(wzl_double,wxl_double);

for ix=1:wxl_double
XW(:,ix)=(ix-wxl_single-1);
end
for iz=1:wzl_double
ZW(iz,:)=(iz-wzl_single-1);
end


for iz=wzl_single+1:nz+wzl_single
    for ix=wxl_single+1:nx+wxl_single
        a=new_img(iz-wzl_single:iz+wzl_single,ix-wxl_single:ix+wxl_single);
        
        A=sum(a(:));
        mean_x=sum(sum(a.*XW))/(A+eps);
        mean_z=sum(sum(a.*ZW))/(A+eps);
        cxx=sum(sum(sum((XW-mean_x).*(XW-mean_x).*a)))/(A+eps);
        czz=sum(sum(sum((ZW-mean_z).*(ZW-mean_z).*a)))/(A+eps);
        cxz=sum(sum(sum((XW-mean_x).*(ZW-mean_z).*a)))/(A+eps);
        CM=[cxx,cxz;cxz,czz];
        [V,D]=eig(CM);

        distt(iz,ix)=abs(V(1,1)*mean_x+V(2,1)*mean_z);
        D1(iz,ix)=D(2,2);
        D2(iz,ix)=D(1,1);
        dipp(iz,ix)=atan2(V(2,2),V(1,2))*180/pi;

    end
end
nmd1=D1(wzl_single+1:nz+wzl_single,wxl_single+1:nx+wxl_single);
nmd2=D2(wzl_single+1:nz+wzl_single,wxl_single+1:nx+wxl_single);
dip=dipp(wzl_single+1:nz+wzl_single,wxl_single+1:nx+wxl_single);
dist=distt(wzl_single+1:nz+wzl_single,wxl_single+1:nx+wxl_single);


