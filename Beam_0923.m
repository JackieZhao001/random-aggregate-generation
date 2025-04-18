clear;
clc;
% % �ߴ�ΪA0=400,B0=100,h0=100�� �������Լ�,
% % ���ǽ����(Ԥ��)
% % 
% % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % !!!!!!!!!!���Ǵ����!!!!!!!!!!!!!!
% % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
num_p=17;
num_2=500;
O2=zeros(num_2,3);         %1 by 3 zeros vector
R2=zeros(num_2,1);
Po2=zeros(num_p,3,num_2);
Fun2=zeros(num_2,30,5);
k=1;
Rat2=0;
Vsum=0;
ScaR=0.10;          %������ȱ�
a0=400;b0=100;h0=100;        %Բ����ߴ�
temp_a=a0-4;
temp_b=b0-4;
temp_h=h0-4;
Vcub=a0*b0*h0;
options=optimset('MaxFunEvals',200000000);
%
% ����������Χ
t_r = 1;
Rat = zeros(1,9);
Pc_d_5 = Fuller(5,20);
Pc_d_6 = Fuller(6,20);
Pc_d_8 = Fuller(8,20);
Pc_d_10 = Fuller(10,20);
Pc_d_12 = Fuller(12,20);
Pc_d_14 = Fuller(14,20);
Pc_d_16 = Fuller(16,20);
Pc_d_18 = Fuller(18,20);
Pc_d_20 = Fuller(20,20);
d_list= [5,6,8,10,12,14,16,18,20];
tic;                            %tic;;;;;;----------------
while Rat2<=0.3
    flag=0;
    m=0;
    psus=zeros(10,1);
    O2(k,1)=unifrnd(-temp_a,temp_a);
    O2(k,2)=unifrnd(-temp_b,temp_b);
    O2(k,3)=unifrnd(4,temp_h);%��������
    if Rat(t_r)< Pc_d_20-Pc_d_18
        R2(t_r)=unifrnd(9,10);%�뾶
    elseif (Rat(t_r)>Pc_d_20-Pc_d_18) && (Rat(t_r)<=Pc_d_20-Pc_d_16)
        R2(t_r)=unifrnd(8,9);
    elseif (Rat(t_r)>Pc_d_20-Pc_d_16) && (Rat(t_r)<=Pc_d_20-Pc_d_14)
        R2(t_r)=unifrnd(7,8);
    elseif (Rat(t_r)>Pc_d_20-Pc_d_14) && (Rat(t_r)<=Pc_d_20-Pc_d_12)
        R2(t_r)=unifrnd(6,7);
    elseif (Rat(t_r)>Pc_d_20-Pc_d_12) && (Rat(t_r)<=Pc_d_20-Pc_d_10)
        R2(t_r)=unifrnd(5,6);
    elseif (Rat(t_r)>Pc_d_20-Pc_d_10) && (Rat(t_r)<=Pc_d_20-Pc_d_8)
        R2(t_r)=unifrnd(4,5);
    else
        break;
    end
%     
%     R2(k)=(1+ScaR)*R2(k);           %(��������Ϊ����)
    %�ж������Ƿ��ڱ߽���
    if O2(k,3)<R2(k) || O2(k,3)>h0-R2(k) || abs(O2(k,1))+R2(k)>a0/2 || abs(O2(k,2))+R2(k)>b0/2
        continue;                             %��������ѭ��
    end   
    %��С���ɷ�Χ
    for j=1:(k-1)
        P=sqrt((O2(k,1)-O2(j,1))^2+(O2(k,2)-O2(j,2))^2+(O2(k,3)-O2(j,3))^2);
        p1=15/20*(R2(k)+R2(j));
        p2=R2(k)+R2(j);
        if P>=p1
            if P<=p2
                m=m+1;
                psus(m)=j;
            end
        else
            flag=1;
            break;
        end
    end
    if flag==1                %��������ѭ��
        continue;
    end
    %     ����17��������
    Po2(:,:,k)=Kpt4(O2(k,1),O2(k,2),O2(k,3),R2(k));
    
    %�ȴ��ж�
    ppsus=zeros(20,1);
    s=0;
    for j=1:m
        n=psus(j);
        %�ж϶������ཻ-----------��һ���б�(k�Ķ����Ƿ���n��)
        Ftemp=reshape(Fun2(n,:,:),30,5);
        Ftemp=Ftemp(:,1:4);
        Ia=Intersect_1(O2(n,:),R2(n),Ftemp,Po2(:,:,k));
        if Ia==1
            flag=1;
            break;
        end
        s=s+1;
        ppsus(s)=n;
    end
    if flag==1                %�ཻ����������ѭ��
        continue;
    end
        %���淽��ϵ��
        Del=delaunayTriangulation(Po2(:,:,k));
        [F0,P0]=freeBoundary(Del);
    for i=1:30
        v1=P0(F0(i,2),:)-P0(F0(i,1),:);           %iΪk�ſ��� ��i����
        v2=P0(F0(i,3),:)-P0(F0(i,1),:);           %1��2��3Ϊ����1��2��3�ŵ�
        v0=cross(v1,v2);
        Fun2(k,i,1)=v0(1); Fun2(k,i,2)=v0(2); Fun2(k,i,3)=v0(3);          %1��2��3\4\5Ϊ����ϵ��
        Fun2(k,i,4)=-(v0(1)*P0(F0(i,1),1)+v0(2)*P0(F0(i,1),2)+v0(3)*P0(F0(i,1),3));
%         Fun2(k,i,5)=(1-ScaR)*Fun2(k,i,4) - ScaR* (v0(1)*O2(k,1) + v0(2)*O2(k,2) + v0(3)*O2(k,3));        
    end
         %�ж϶������ཻ-----------�ڶ����б�(n�Ķ����Ƿ���k��)
    for j=1:s
        n=ppsus(j);
        Ftemp=reshape(Fun2(k,:,:),30,5);
        Ftemp=Ftemp(:,1:4);
        Ia=Intersect_1(O2(k,:),R2(k),Ftemp,Po2(:,:,n));
        if Ia==1
            flag=1;
            break;
        end
        %     �ж϶������ཻ-----------�������б�(�Ƿ��ཻ)
        Ia=intersectionHull4(1,Po2(:,:,k),1,Po2(:,:,n),options);
        if Ia==1
            flag=1;
            break;
        end
    end
    if flag==1                %�ཻ����������ѭ��
        continue;
    end
    %���������
    [~,v]=convexHull(Del);
    Vsum=Vsum+v;
    Rat2=Vsum/Vcub;
    t_r = t_r+1;
    Rat(t_r)=Rat2;
    
    %     %����������ĵ㵽���ľ���
    %     sumR=0;
    %     for i=1:12
    %         cen=1/3*(P0(F0(i,3),:)+P0(F0(i,2),:)+P0(F0(i,1),:));           %��i���� ������
    %         d=sqrt((cen(1)-O2(k,1))^2+(cen(2)-O2(k,2))^2+(cen(3)-O2(k,3))^2) ;
    %         sumR=sumR+d/R2(k);
    %     end
    %     avr(k)=sumR/12;
    
    fprintf('%d***********%d \n',k,Rat2);
    k=k+1;
end
% toc;
% save('cylbd_91_1.mat');


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!������С����!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

num_1=2000;
Fun1=zeros(num_1,30,5);
O1=zeros(num_1,3);
R1=zeros(num_1,1);
Po1=zeros(num_p,3,num_1);
kk=1;
Vsum=0;
Rat1=0;
ScaR=0.2;
temp_a=a0-2.5;
temp_b=b0-2.5;
temp_h=h0-2.5;
tic;                            %tic;;;;;;----------------
while Rat1<=0.065
    flag=0;
    m=0;
    mm=0;
    psus=zeros(20,2);
    O1(kk,1)=unifrnd(-temp_a,temp_a);
    O1(kk,2)=unifrnd(-temp_b,temp_b);
    O1(kk,3)=unifrnd(2.5,temp_h);%��������
    
    if  Rat(t_r)<=Pc_d_20-Pc_d_6
        R1(kk)=unifrnd(3,4);
    elseif (Rat(t_r)>Pc_d_20-Pc_d_6) && (Rat(t_r)<=Pc_d_20-Pc_d_5)
        R1(kk)=unifrnd(2.5,3);
    else
        break;
    end
    
    %�ж������Ƿ��ڱ߽���
    if O1(kk,3)<R1(kk) || O1(kk,3)>h0-R1(kk) || abs(O1(kk,1))+R1(kk)>a0/2 || abs(O1(kk,2))+R1(kk)>b0/2
        continue;                             %��������ѭ��
    end
    %��С�����ɷ�Χ------Po2
    for j=1:k-1
        P=sqrt((O1(kk,1)-O2(j,1))^2+(O1(kk,2)-O2(j,2))^2+(O1(kk,3)-O2(j,3))^2);
        p1=13/20*R2(j)+R1(kk);
        p2=R1(kk)+R2(j);
        if P>=p1
            if P<=p2
                mm=mm+1;
                psus(mm,2)=j;
            end
        else
            flag=1;
            break;
        end
    end
    if flag==1                 %��������ѭ��
        continue;
    end
    %��С�����ɷ�Χ-------Po1
    for j=1:(kk-1)
        P=sqrt((O1(kk,1)-O1(j,1))^2+(O1(kk,2)-O1(j,2))^2+(O1(kk,3)-O1(j,3))^2);
        p1=13/20*R1(j)+R1(kk);
        p2=R1(kk)+R1(j);
        if P>=p1
            if P<=p2
                m=m+1;
                psus(m,1)=j;
            end
        else
            flag=1;
            break;
        end
    end
    if flag==1               %��������ѭ��
        continue;
    end
    %     ����12��������
    Po1(:,:,kk)=Kpt4(O1(kk,1),O1(kk,2),O1(kk,3),R1(kk));
    %�ȴ��ж�
    %�ж϶������ཻ-----------��Po2��һ���б�(k�Ķ����Ƿ���n��)
    for j=1:mm
        n=psus(j,2);
        Ftemp=reshape(Fun2(n,:,:),30,5);
        Ftemp=Ftemp(:,1:4);
        Ia=Intersect_1(O2(n,:),R2(n),Ftemp,Po1(:,:,kk));
        if Ia==1
            flag=1;
            break;
        end
    end
    if flag==1                %�ཻ����������ѭ��
        continue;
    end
    for j=1:m
        n=psus(j,1);
        %�ж϶������ཻ-----------��Po1��һ���б�(k�Ķ����Ƿ���n��)
        Ftemp=reshape(Fun1(n,:,:),30,5);
        Ftemp=Ftemp(:,1:4);
        Ia=Intersect_1(O1(n,:),R1(n),Ftemp,Po1(:,:,kk));
        if Ia==1
            flag=1;
            break;
        end
    end
    if flag==1                %�ཻ����������ѭ��
        continue;
    end
    %���淽��ϵ��
    Del=delaunayTriangulation(Po1(:,:,kk));
    [F0,P0]=freeBoundary(Del);
    for i=1:30
        v1=P0(F0(i,2),:)-P0(F0(i,1),:);           %iΪk�ſ��� ��i����
        v2=P0(F0(i,3),:)-P0(F0(i,1),:);           %1��2��3Ϊ����1��2��3�ŵ�
        v0=cross(v1,v2);
        Fun1(kk,i,1)=v0(1); Fun1(kk,i,2)=v0(2); Fun1(kk,i,3)=v0(3);          %1��2��3\4\5Ϊ����ϵ��
        Fun1(kk,i,4)=-(v0(1)*P0(F0(i,1),1)+v0(2)*P0(F0(i,1),2)+v0(3)*P0(F0(i,1),3));
%       Fun1(kk,i,5)=-ScaR*Fun1(kk,i,4) - ScaR*(v0(1)*O1(kk,1)+v0(2)*O1(kk,2)+v0(3)*O1(kk,3));
    end
    %�ж϶������ཻ-----------Po2�ڶ����б�(n�Ķ����Ƿ���kk��)
    for j=1:mm
        n=psus(j,2);
        Ftemp=reshape(Fun1(kk,:,:),30,5);
        Ftemp=Ftemp(:,1:4);
        Ia=Intersect_1(O1(kk,:),R1(kk),Ftemp,Po2(:,:,n));
        if Ia==1
            flag=1;
            break;
        end
    end
    if flag==1                %�ཻ����������ѭ��
        continue;
    end
    %�ж϶������ཻ-----------Po1�ڶ����б�(n�Ķ����Ƿ���kk��)
    for j=1:m
        n=psus(j,1);
        Ftemp=reshape(Fun1(kk,:,:),30,5);
        Ftemp=Ftemp(:,1:4);
        Ia=Intersect_1(O1(kk,:),R1(kk),Ftemp,Po1(:,:,n));
        if Ia==1
            flag=1;
            break;
        end
    end
    if flag==1                %�ཻ����������ѭ��
        continue;
    end
    %     �ж϶������ཻ-----------Po2�������б�(�Ƿ��ཻ)
    for j=1:mm
        n=psus(j,2);
        Ia=intersectionHull4(1,Po1(:,:,kk),1,Po2(:,:,n),options);
        if Ia==1
            flag=1;
%             fprintf('ok22222222222ok\n');
            break;
        end
    end
    if flag==1                %�ཻ����������ѭ��
        continue;
    end
        %     �ж϶������ཻ-----------Po1�������б�(�Ƿ��ཻ)
    for j=1:m
        n=psus(j,1);
        Ia=intersectionHull4(1,Po1(:,:,kk),1,Po1(:,:,n),options);
        if Ia==1
            flag=1;
            fprintf('ok11111111111ok\n');
            break;
        end
    end
    if flag==1                %�ཻ����������ѭ��
        continue;
    end 
    %���������
    [~,v]=convexHull(Del);
    Vsum=Vsum+v;
    Rat1=Vsum/Vcub;
    t_r = t_r+1;
    Rat(t_r)=Rat1;
    fprintf('%d******************************* %d \n',kk,Rat1);
    kk=kk+1;
end
% toc;
save('Beam_0928_2.mat')
