ug1=0;
ug2=0;
ug3=0;
ug0=0;
ug4=0;
yg0=[0;0.2;0.6];
yg1=[0;0.2;0.6];
yg2=yg1;
yg3=yg2;
yg4=yg3;
r=[0.4;0.3];
Rd=0.001;
Qd=[0.001 0 0;
    0 0.1 -0.1;
    0 -0.1 0.1];
gamma=0.575;
siold=[0;0.2;0.6;0.6];
Pd_=ones(24)*0.0001;
rpre=[r;r;r;r];
rnext=rpre;
ukdata=zeros(5000,1);
ykdata=zeros(3,5000);
test=[];
for k=1:5000
    upre1=[ug1;ug2;ug3];
    ypre=[yg0;yg1;yg2;yg3];
    ek=nhieu(k);
    ug0=-1/(Rd/gamma+Pd_(1,1))*(Pd_(1,2:4)*upre1+Pd_(1,5:16)*ypre+Pd_(1,17:24)*rpre)+ek;
    yg4=yg3;
    yg3=yg2;
    yg2=yg1;
    yg1=yg0;
    [yg0,sinew]=mohinh(siold,ug0,r);
    ug4=ug3;
    ug3=ug2;
    ug2=ug1;
    ug1=ug0;
    siold=sinew;
    ukdata(k)=ug0;
    ykdata(:,k)=yg0;
end
Odata=zeros(4996,576);
for k=5:5000
    upre=ukdata(k-1:-1:k-4);
    yprestk=ykdata(:,k-1:-1:k-4);
    ypre=reshape(yprestk,[12 1]);
    zpre=[upre;ypre;rpre];
    Odata(k-4,:)=kron(zpre',zpre');
end
ujdata=ukdata;
for j=0:1000
    Vdata=zeros(4996,1);
    for k=5:5000
        unext=ukdata(k:-1:k-3);
        ynextstk=ykdata(:,k:-1:k-3);
        ynext=reshape(ynextstk,[12 1]);
        znext=[unext;ynext;rnext];
        Vdata(k-4)=ykdata(:,k)'*Qd*ykdata(:,k)+ujdata(k)'*Rd*ujdata(k)+gamma*znext'*Pd_*znext;
    end
    Pnewstk=pinv(Odata)*Vdata;
    Pnew=reshape(Pnewstk,[24 24]);
    if (norm(Pnew-Pd_,1)<0.0001)
        break;
    end
    for k=5:5000
        upre1=ukdata(k-1:-1:k-3);
        ynextstk=ykdata(:,k:-1:k-3);
        ynext=reshape(ynextstk,[12 1]);
        ujdata(k)=-1/(Rd/gamma+Pnew(1,1))*(Pnew(1,2:4)*upre1+Pnew(1,5:16)*ynext+Pnew(1,17:24)*rnext);
    end
    Pd_=Pnew;
    test=[test Vdata];
end