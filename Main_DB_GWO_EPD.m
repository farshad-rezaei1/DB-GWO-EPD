%_________________________________________________________________________%
% DB-GWO-EPD: A Grey Wolf Optimizer Equipped with Diversity-Based         %
% Evolutionary Population Dynamics                                        %
%                                                                         %
% Developed in MATLAB R2018b                                              %
%                                                                         %
% Inventor and programmer: Farshad Rezaei, PhD                            %
%                                                                         %
% e-Mail: farshad.rezaei@gmail.com                                        %
%         f.rezaei@alumni.iut.ac.ir                                       %
%                                                                         %
% Homepage: https://www.linkedin.com/in/farshad-rezaei-5a92559a/          %
%                                                                         %
% Main paper: Rezaei, F.; Safavi, H.R.; Abd Elaziz, M.; Abualigah, L.;    %
% Mirjalili, S.; Gandomi, A.H. Diversity-Based Evolutionary Population    %
% Dynamics: A New Operator for Grey Wolf Optimizer. Processes 2022, 10,   %
% 2615. https://doi.org/10.3390/pr10122615                                %
%_________________________________________________________________________%

% The initial parameters that you need are:
%_________________________________________________________________________
% fobj=@YourCostFunction
% nx=number of your variables
% lb=the lower bound of variables which can generally be a fixed number or a vector
% ub=the upper bound of variables which can generally be a fixed number or a vector
% notice: if the lower nad upper bounds are not fixed for all variables, 
% they appear in the forms of the vectors "varmin" and "varmax", as illustrated in following

% To run DB-GWO-EPD: [z_iter,z_final,pos_final]=DB_GWO_EPD(np,nx,maxit,varmax,varmin,a_max,a_min,mut_max,mut_min,fobj)

%_________________________________________________________________________
% Set the required parameters to run the DB-GWO-EPD algorithm

% This code is for solving the minimization problems. To maximize a desired 
% cost function,please implement this code upon inverting the sign of the cost function

% function [z_best]=Main_DB_GWO_EPD(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num) % Only for CEC2017
clc
clear
close all
tic
Function_name='F2'; % Name of the test function that can be from 'F1' to 'F13'
Function_set=1; % 1: Standard, 2:  CEC2017
if Function_set==1
    np=30; % Number of solutions
    maxit=1000; % Maximum number of iterations
    [lb,ub,nx,fobj]=Objective_Function_DB_GWO_EPD(Function_name); % Load details of the selected benchmark function
elseif Function_set==2
    np=pop_size;
    maxit=iter_max;
    lb=Xmin;
    ub=Xmax;
    nx=D;
end 
run=30; % Maximum number of the algorithm runnings conducted
a_max=2; % Upper bound of the acceleration coefficient
a_min=0; % Lower bound of the acceleration coefficient
mut_max=1; % Upper bound of the mutation step
mut_min=0; % Lower bound of the mutation step
varmax=ub*ones(1,nx); % Upper bound defined for the positions which can generally be a desired vector
varmin=lb*ones(1,nx); % Lower bound defined for the positions which can generally be a desired vector
z_iter_main=zeros(run,maxit);
z_final_main=zeros(run);
pos_final_main=zeros(run,nx);
x1=zeros(maxit);
y1=zeros(maxit);

% Run the DB-GWO-EPD algorithm for "run" times 
for nrun=1:run
    [z_iter,z_final,pos_final]=DB_GWO_EPD(np,nx,maxit,varmax,varmin,a_max,a_min,mut_max,mut_min,fobj,Function_set); % for Standard Functions
%     [z_iter,z_final,pos_final]=DB_GWO_EPD(np,nx,maxit,varmax,varmin,a_max,a_min,mut_max,mut_min,fhd,func_num,Function_set); % for CEC2017
    z_iter_main(nrun,1:maxit)=z_iter(1:maxit);
    z_final_main(nrun)=z_final;
    pos_final_main(nrun,1:nx)=pos_final(1:nx);
end
z_best=z_final;

% Display the comprehensive results
if Function_set==1
    disp(['The final statistical results calculated when implementing the DB-GWO-EPD algorithm for ',num2str(run),' times are as follows:']);
    disp(['The average of the final objective function values calculated over ',num2str(run),' times = ',num2str(mean(z_final_main(1:run)))]);
    % disp(['The median of the final objective function values calculated over ',num2str(run),' times = ',num2str(median(z_final_main(1:run)))]);
    % disp(['The best of the final objective function values calculated over ',num2str(run),' times = ',num2str(min(z_final_main(1:run)))]);
    disp(['The standard deviation of the final objective function values calculated over ',num2str(run),' times = ',num2str(std(z_final_main(1:run)))]);
end

% Plot the convergence curve of the DB-GWO-EPD over the course of iterations
for i=1:maxit
    x1(i)=i;sum1=0;
    for j=1:run
        sum1=sum1+z_iter_main(j,i);
    end
    y1(i)=sum1/run;
end
semilogy(x1,y1,'-r')
xlabel('Iteration');
ylabel('Average best-so-far');
legend('DB-GWO-EPD');
hold on
time_db_gwo_epd = toc;
disp(['Elapsed time of running the DB-GWO-EPD for ',num2str(run),' times = ',num2str(time_db_gwo_epd),' seconds']);