% clear all
% mex cec17_func.cpp -DWINDOWS
% func_num=1;
% clc
% clear
% close all
tic
flag=1;
D=50;
Xmin=-100;
Xmax=100;
pop_size=50;
iter_max=1000;
fes=iter_max*pop_size;
runs=30;
z=zeros(runs);
runtime=zeros(runs);
fhd=str2func('cec17_func');
for i=1:1
func_num=i;
if func_num == 11 || func_num == 21
    disp(' ');
end
for j=1:runs
% DB-GWO-EPD
if flag==1
    [z_best]=Main_DB_GWO_EPD(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
    z(j)=z_best;
end
end
disp(num2str(mean(z(1:runs))));
% disp(num2str(median(z(1:runs))));
% disp(num2str(min(z(1:runs))));
disp(num2str(std(z(1:runs))));
end
toc