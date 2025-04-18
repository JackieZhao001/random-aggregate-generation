clear;
clc;
%处理生成不带界面区的 骨料单元编号
load 'cyl_241104_1.mat'
% 
O2(k:num_2,:)=[];
O1(kk:num_1,:)=[];

R2(k:num_2)=[];
R1(kk:num_1)=[];
O=[O2;O1];
R=[R2;R1];
% 
Fun2(k:num_2,:,:)=[];
Fun1(kk:num_1,:,:)=[];
num=k+kk-2;
Fun3=zeros(num,30,5);
Fun3(1:k-1,:,:)=Fun2;
Fun3(k:num,:,:)=Fun1;
% 
Po2(:,:,k:num_2)=[];
Po1(:,:,kk:num_1)=[];
Po3=zeros(17,3,num);
Po3(:,:,1:k-1)=Po2;
Po3(:,:,k:num)=Po1;
% 
% Fun3=Add_Bound(Po3,k,Fun3,O);

eloc=importdata('Elemloc_20_4.txt');
eloc=eloc.*1000;
ecount=length(eloc);
l0=1/2*ecount;
emat2=zeros(l0,1);
emat3=zeros(l0,1);
enum=zeros(ecount,1); %单元编号
for ii=1:ecount
    enum(ii)=ii;
end
ecount0=ecount;
tic;
l2=0;    %elements of mat=2  number（界面层）
l3=0;   %elements of mat=3  number(骨料)
for ii=1:num
    S=R(ii)^2; %球形范围内
    ecount=ecount0;
    for jj=1:ecount
        D=(eloc(jj,1)-O(ii,1))^2+(eloc(jj,2)-O(ii,2))^2+(eloc(jj,3)-O(ii,3))^2;
        if D>S
            continue;       %判断下一个单元
        end
        flag=0;
        fcheck=0;           %仅需要判断out时，fcheck为1；需要判断in时，fcheck为0。
        for qq=1:30
            F_in=Fun3(ii,qq,1)*eloc(enum(jj),1)+Fun3(ii,qq,2)*eloc(enum(jj),2)+Fun3(ii,qq,3)*eloc(enum(jj),3)+Fun3(ii,qq,4);
            F_out=F_in-Fun3(ii,qq,5);
            if F_out>0
                flag=1;
                break;
            elseif fcheck==0
                if F_in>0
                   fcheck=1;
                end
            end
        end
        if flag==1
            continue;
        end
        if fcheck==1
            l2=l2+1;
            emat2(l2)=jj;     %yes,note the element
        else
            l3=l3+1;
            emat3(l3)=jj;     %yes,note the element
        end
      
        fprintf("i=%d \n",ii);
%         
%         temp1=eloc(enum(jj),:); %exchange the eloc
%         eloc(enum(jj),:)=eloc(enum(ecount0),:);
%         eloc(enum(ecount0),:)=temp1;
%         
%         temp2=enum(jj);         %exchange the enum
%         enum(jj)=enum(ecount);
%         enum(ecount)=temp2;
%         
%         ecount0=ecount0-1;      %so the next aggergate need not judge this elements.
    end
end
toc;       


l3=l3+1;
emat3(l3:l0)=[];
emat3=unique(emat3);
[l3,~]=size(emat3);
% %
% %
Ri=intersect(emat2(:),emat3(:));              % 求两向bai量（矩阵）的交du集
[l,~]=size(Ri);
emat32=zeros(l2-l,1);
j=1;
for i=1:l3                                               %骨料
    s=emat3(i);
    a=ismember(s,Ri);
    if a==0
        emat32(j)=s;
        j=j+1;
    end
     
end

emat3=emat32';
% 
EN=jj-l3;
emat0=zeros(EN,1);
j=1;
for i=1:jj
    bb=ismember(i,emat3);
%     if double(bb)==0
%         emat0(j)=i;
%         j=j+1;
%     end
    if double(bb)==1
        continue;
    else
        emat0(j)=i;
        j=j+1;
    end
    if bb==1
        fprintf('wrong=%d\n','i');
    end
end

% 
fid = fopen('D:\workfile\Matlab\Guliao_new\D0_ematA_1104.txt','wt');
fprintf(fid,'%d\n',emat3);
fclose(fid);
fid = fopen('D:\workfile\Matlab\Guliao_new\D0_ematC_1104.txt','wt');
fprintf(fid,'%d\n',emat0);





