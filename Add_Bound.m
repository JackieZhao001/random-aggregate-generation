function Fun3=Add_Bound(Po3,k,Fun3,O)
% 添加界面层——添加fun（：，：，5）的值
[~,~,num]=size(Po3);
for j=1:1:k-1
    ScaR=0.1;
    Del=delaunayTriangulation(Po3(:,:,j));
    [F0,P0]=freeBoundary(Del);
    for i=1:30
        v1=P0(F0(i,2),:)-P0(F0(i,1),:);           %i为k号颗粒 第i个面
        v2=P0(F0(i,3),:)-P0(F0(i,1),:);           %1、2、3为面上1，2，3号点
        v0=cross(v1,v2);
%         Fun3(j,i,1)=v0(1); Fun3(j,i,2)=v0(2); Fun3(j,i,3)=v0(3);          %1、2、3\4\5为面上系数
%         Fun3(j,i,4)=-(v0(1)*P0(F0(i,1),1)+v0(2)*P0(F0(i,1),2)+v0(3)*P0(F0(i,1),3));
        Fun3(j,i,5)=-ScaR*Fun3(j,i,4) - ScaR*(v0(1)*O(j,1)+v0(2)*O(j,2)+v0(3)*O(j,3));
    end
end
ScaR=0.25;
for j=k-1:1:num
    Del=delaunayTriangulation(Po3(:,:,j));
    [F0,P0]=freeBoundary(Del);
    for i=1:30
        v1=P0(F0(i,2),:)-P0(F0(i,1),:);           %i为k号颗粒 第i个面
        v2=P0(F0(i,3),:)-P0(F0(i,1),:);           %1、2、3为面上1，2，3号点
        v0=cross(v1,v2);
%         Fun3(j,i,1)=v0(1); Fun3(j,i,2)=v0(2); Fun3(j,i,3)=v0(3);          %1、2、3\4\5为面上系数
%         Fun3(j,i,4)=-(v0(1)*P0(F0(i,1),1)+v0(2)*P0(F0(i,1),2)+v0(3)*P0(F0(i,1),3));
        Fun3(j,i,5)=-ScaR*Fun3(j,i,4) - ScaR*(v0(1)*O(j,1)+v0(2)*O(j,2)+v0(3)*O(j,3));
    end
end
