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

% This function is to initialize the position and velocity of the wolves to start the optimization process
function [pp]=Initialization_DB_GWO_EPD(np,nx,varmax,varmin)
pp=zeros(np,nx); 
for j=1:np
    pp(j,1:nx)=(varmax-varmin).*rand(1,nx)+varmin;
end