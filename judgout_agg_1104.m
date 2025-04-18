clear;  % 清除工作区
clc;    % 清除命令行窗口

% 处理生成不带界面区的骨料单元编号
% load 'cyl_241104_1.mat';  % 加载数据文件
load 'cyl0516_1.mat';  % 加载数据文件cyl0516_1

% 删除不需要的部分数据
% O2(k:num_2,:) = [];  % 删除O2矩阵中的行
% O1(kk:num_1,:) = [];  % 删除O1矩阵中的行
% 
% R2(k:num_2) = [];  % 删除R2数组中的元素
% R1(kk:num_1) = [];  % 删除R1数组中的元素

% 合并处理后的O和R
% O = [O2; O1];  % 合并O2和O1
% R = [R2; R1];  % 合并R2和R1

% 删除不需要的Fun数据
% Fun2(k:num_2,:,:) = [];  % 删除Fun2中的行
% Fun1(kk:num_1,:,:) = [];  % 删除Fun1中的行

% num = k + kk - 2;  % 计算有效单元数量
% Fun3 = zeros(num, 30, 5);  % 初始化Fun3矩阵
% Fun3(1:k-1,:,:) = Fun2;  % 将Fun2的值赋给Fun3
% Fun3(k:num,:,:) = Fun1;  % 将Fun1的值赋给Fun3
% 
% % 删除不需要的Po数据
% Po2(:,:,k:num_2) = [];  % 删除Po2中的列
% Po1(:,:,kk:num_1) = [];  % 删除Po1中的列
% 
% % 合并处理后的Po数据
% Po3 = zeros(17, 3, num);  % 初始化Po3矩阵
% Po3(:,:,1:k-1) = Po2;  % 将Po2的值赋给Po3
% Po3(:,:,k:num) = Po1;  % 将Po1的值赋给Po3

% Fun3 = Add_Bound(Po3, k, Fun3, O);  % 调用Add_Bound函数（已注释）

% 导入单元位置数据
eloc = importdata('Elemloc_20_4_h20.txt');  % 导入数据
eloc = eloc .* 1000;  % 将位置数据缩放
ecount = length(eloc);  % 获取单元数量
l0 = 1/2 * ecount;  % 计算l0
emat2 = zeros(l0, 1);  % 初始化emat2
emat3 = zeros(l0, 1);  % 初始化emat3
enum = zeros(ecount, 1);  % 初始化单元编号数组

% 为enum数组赋值
for ii = 1:ecount
    enum(ii) = ii;  % 赋值为1到ecount
end

ecount0 = ecount;  % 备份ecount
tic;  % 开始计时

l2 = 0;  % 材料2的单元数量（界面层）
l3 = 0;  % 材料3的单元数量（骨料）

% 遍历每个单元
for ii = 1:num
    S = R(ii)^2;  % 球形范围内的半径平方
    ecount = ecount0;  % 恢复ecount值

    % 遍历每个位置
    for jj = 1:ecount
        D = (eloc(jj,1) - O(ii,1))^2 + (eloc(jj,2) - O(ii,2))^2 + (eloc(jj,3) - O(ii,3))^2;  % 计算距离
        if D > S
            continue;  % 如果距离超出范围，跳过
        end
        
        flag = 0;  % 标记变量
        fcheck = 0;  % 检查变量

        % 遍历Fun3中的30个特征
        for qq = 1:30
            F_in = Fun3(ii,qq,1) * eloc(enum(jj),1) + Fun3(ii,qq,2) * eloc(enum(jj),2) + Fun3(ii,qq,3) * eloc(enum(jj),3) + Fun3(ii,qq,4);  % 计算输入
            F_out = F_in - Fun3(ii,qq,5);  % 计算输出

            if F_out > 0
                flag = 1;  % 如果输出大于0，设置标记
                break;  % 退出循环
            elseif fcheck == 0
                if F_in > 0
                   fcheck = 1;  % 如果输入大于0，设置检查标记
                end
            end
        end
        
        if flag == 1
            continue;  % 如果标记为1，跳过
        end

        if fcheck == 1
            l2 = l2 + 1;  % 增加材料2的计数
            emat2(l2) = jj;  % 记录该元素
        else
            l3 = l3 + 1;  % 增加材料3的计数
            emat3(l3) = jj;  % 记录该元素
        end
      
        fprintf("i=%d \n", ii);  % 打印当前索引
    end
end
toc;  % 结束计时

l3 = l3 + 1;  % 更新l3
emat3(l3:l0) = [];  % 删除emat3中多余的部分
emat3 = unique(emat3);  % 取emat3中的唯一值
[l3,~] = size(emat3);  % 获取emat3的大小

% 计算交集
Ri = intersect(emat2(:), emat3(:));  % 求emat2和emat3的交集
[l,~] = size(Ri);  % 获取交集的大小
emat32 = zeros(l2-l, 1);  % 初始化emat32
j = 1;

% 遍历骨料单元
for i = 1:l3
    s = emat3(i);
    a = ismember(s, Ri);  % 检查是否在交集中
    if a == 0
        emat32(j) = s;  % 如果不在，记录该元素
        j = j + 1;
    end
end

emat3 = emat32';  % 转置emat3

% 计算不在emat3中的元素数量
EN = jj - l3;  
emat0 = zeros(EN, 1);  % 初始化emat0
j = 1;

% 遍历所有元素
for i = 1:jj
    bb = ismember(i, emat3);  % 检查是否在emat3中
    if double(bb) == 1
        continue;  % 如果在，跳过
    else
        emat0(j) = i;  % 如果不在，记录该元素
        j = j + 1;
    end
    if bb == 1
        fprintf('wrong=%d\n','i');  % 打印错误信息
    end
end

% 保存结果到文本文件
% 新建文件并写入emat3
fid = fopen('E:\Work_file\Matlab\1_Mine\Guliao_new\D0_ematA_1104.txt', 'w'); % 使用 'w' 模式
if fid == -1
    error('无法打开文件 D0_ematA_1104.txt，请检查路径和权限。');
end
fprintf(fid, '%d\n', emat3);
fclose(fid);

% 新建文件并写入emat0
fid = fopen('E:\Work_file\Matlab\1_Mine\Guliao_new\D0_ematC_1104.txt', 'w'); % 使用 'w' 模式
if fid == -1
    error('无法打开文件 D0_ematC_1104.txt，请检查路径和权限。');
end
fprintf(fid, '%d\n', emat0);
fclose(fid);

