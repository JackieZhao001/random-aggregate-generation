% clear;
% clc;
% 
% emat1count = length(emat0);
% emat2count = length(emat2);
% emat3count = length(emat32);
% full
% emat1count=11934;
% emat2count=2975;
% emat3count=9091;
emat1count=95838;
emat2count=23865;
emat3count=72297;

NEMAX=emat1count+emat3count+emat2count;
% %%%%%%%%%%%%弹性模量
% % 
%    砂浆
format long;
digits(6);
B=3;  %%%top
EC=37e9; %
ER1=roundn(wblrnd(EC,B,[emat1count,1]),4);
%    ITZ
B=3;%%%coal
Eitz=29e9; %
ER2=roundn(wblrnd(Eitz,B,[emat2count,1]),4);
%    骨料
B=4;%%%coal
EA=70e9; %
ER3=roundn(wblrnd(EA,B,[emat3count,1]),4);
%
Emean=mean(ER1)*0.5+mean(ER3)*0.5;
%
% 
fid = fopen('E:\Work_file\Matlab\1_Mine\Guliao_new\D1_WblEx1.txt','wt');
fprintf(fid,'%g\n',ER1);
fclose(fid);
%
fid = fopen('E:\Work_file\Matlab\1_Mine\Guliao_new\D1_WblEx2.txt','wt');
fprintf(fid,'%g\n',ER2);
fclose(fid);
% %
fid = fopen('E:\Work_file\Matlab\1_Mine\Guliao_new\D1_WblEx3.txt','wt');
fprintf(fid,'%g\n',ER3);
fclose(fid);

