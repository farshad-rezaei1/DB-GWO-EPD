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

% DB-GWO-EPD algorithm  
function [z_iter,z_final,pos_final]=DB_GWO_EPD(np,nx,maxit,varmax,varmin,a_max,a_min,mut_max,mut_min,fobj,Function_set)
% function [z_iter,z_final,pos_final]=DB_GWO_EPD(np,nx,maxit,varmax,varmin,a_max,a_min,mut_max,mut_min,fhd,func_num,Function_set) % Only for CEC2017
it=1;
% disp(['Number of Iterations = ',num2str(it)]);
optimal_pos=zeros(maxit,nx);
z=zeros(np);
div=zeros(np);
index=zeros(np);
ss=zeros(np);
pos_final=zeros(nx);
z_iter=zeros(maxit);
z_alpha=inf;
z_beta=inf;
z_delta=inf;
alpha=zeros(1,nx);
beta=zeros(1,nx);
delta=zeros(1,nx);
alpha_div=zeros(1,nx);
beta_div=zeros(1,nx);
delta_div=zeros(1,nx);
nbest=round(np/2);

% Initialization process of the algorithm
[pp]=Initialization_DB_GWO_EPD(np,nx,varmax,varmin);

% Objective function evaluations 
for j=1:np
    if Function_set==1
        z(j)=fobj(pp(j,1:nx)); % for Standard Functions
    elseif Function_set==2
        z(j)=Objective_Function_DB_GWO_EPD(pp(j,1:nx),fhd,func_num); % for CEC2017
    end
end

% Starting the optimization process
for j=1:np
    index(j)=j;
end

% Determining Alpha, Beta, and Delta wolves in the current run
for j=1:np
    if z(j)<=z_alpha
        z_alpha=z(j);
        alpha(1,1:nx)=pp(index(j),1:nx);
    elseif z(j)>z_alpha && z(j)<=z_beta
        z_beta=z(j);
        beta(1,1:nx)=pp(index(j),1:nx);
    elseif z(j)>z_beta && z(j)<=z_delta
        z_delta=z(j);
        delta(1,1:nx)=pp(index(j),1:nx);
    end
end

% Sorting the objective values and their corresponding solutions 
for j=1:np-1
    for jj=j+1:np
        if z(jj)<z(j)
            c1=z(j);
            z(j)=z(jj);
            z(jj)=c1;
            c2=index(j);
            index(j)=index(jj);
            index(jj)=c2;
        end
    end
end

% Calculating the diversity index for each solution (Eq. (14))
for j=1:np
    k=0;
    for jj=1:np
        if (jj~=j) 
            k=k+1;
            ss(j)=abs(z(j)-z(jj));
            if k==1
                min_ss=ss(j);
            elseif ss(j)<min_ss
                min_ss=ss(j);
            end
        end
    end
    if k==0
        div(j)=inf;
    else
        div(j)=min_ss;
    end
end

% Determining three most-diversified solutions in the current run
div_alpha=-inf;
div_beta=-inf;
div_delta=-inf;
for j=1:np
    if div(j)>=div_alpha
        div_alpha=div(j);
        alpha_div(1,1:nx)=pp(index(j),1:nx);
    elseif div(j)<div_alpha && div(j)>=div_beta
        div_beta=div(j);
        beta_div(1,1:nx)=pp(index(j),1:nx);
    elseif div(j)<div_beta && div(j)>=div_delta
        div_delta=div(j);
        delta_div(1,1:nx)=pp(index(j),1:nx);
    end
end

% Repositioning the first half of the best-fitted solutions around the three most diversified ones (Eqs. (15-17))
for j=1:np
    if j<=nbest
        jjj=index(j);
        random=rand(1,1);
        if random<=(1/3)
            guide(1,1:nx)=alpha_div(1,1:nx)+sign(sign(rand(1,1)-0.5)+0.5)*((varmax(1,1:nx)-varmin(1,1:nx)).*rand(1,nx)+varmin(1,1:nx));
        elseif random>(1/3) && random<=(2/3)
            guide(1,1:nx)=beta_div(1,1:nx)+sign(sign(rand(1,1)-0.5)+0.5)*((varmax(1,1:nx)-varmin(1,1:nx)).*rand(1,nx)+varmin(1,1:nx));
        elseif random>(2/3) && random<=1
            guide(1,1:nx)=delta_div(1,1:nx)+sign(sign(rand(1,1)-0.5)+0.5)*((varmax(1,1:nx)-varmin(1,1:nx)).*rand(1,nx)+varmin(1,1:nx));
        end
        pp(jjj,1:nx)=guide;
    end
end
z_optimal(it)=z_alpha;
optimal_pos(it,:)=alpha(1,:);

% Saving the best-so-far objective value in the current run
z_iter(it)=z_optimal(it);

% The Main Loop
while it<maxit
    it=it+1;
%     disp(['Number of Iterations= ',num2str(it)]);
    aa=a_max-(a_max-a_min)*(it-1)/(maxit-1);
    mut=mut_max-(mut_max-mut_min)*(it-1)/(maxit-1);
    for j=1:np 
        a_alpha(1,1:nx)=(2*rand(1,nx)-ones(1,nx))*aa; 
        a_beta(1,1:nx)=(2*rand(1,nx)-ones(1,nx))*aa; 
        a_delta(1,1:nx)=(2*rand(1,nx)-ones(1,nx))*aa; 
        guide_alpha=alpha(1,1:nx)-a_alpha(1,1:nx).*abs(alpha(1,1:nx)-pp(j,1:nx)); % Eq. (6)
        guide_beta=beta(1,1:nx)-a_beta(1,1:nx).*abs(beta(1,1:nx)-pp(j,1:nx)); % Eq. (7)
        guide_delta=delta(1,1:nx)-a_delta(1,1:nx).*abs(delta(1,1:nx)-pp(j,1:nx)); % Eq. (8)
        pp(j,1:nx)=(guide_alpha+guide_beta+guide_delta)/3; % Eq. (19)
            
        % Return back the solution positions if going beyond the position boundaries
        flag4lbp=pp(j,:)<varmin(1,:);
        flag4ubp=pp(j,:)>varmax(1,:);
        varmin_new(1,:)=varmin(1,:)+rand(1,nx).*mut.*(varmax(1,:)-varmin(1,:)); 
        varmax_new(1,:)=varmax(1,:)-rand(1,nx).*mut.*(varmax(1,:)-varmin(1,:));
        pp(j,:)=pp(j,:).*(~(flag4lbp+flag4ubp))+varmin_new.*flag4lbp+varmax_new.*flag4ubp;
    
    % Objective function evaluations 
    if Function_set==1
        z(j)=fobj(pp(j,1:nx)); % for Standard Functions
    elseif Function_set==2
        z(j)=Objective_Function_DB_GWO_EPD(pp(j,1:nx),fhd,func_num); % for CEC2017
    end
    end
    
    % Starting the optimization process
    for j=1:np
        index(j)=j;
    end

    % Determining Alpha, Beta, and Delta wolves in the current run
    for j=1:np
        if z(j)<=z_alpha
            z_alpha=z(j);
            alpha(1,1:nx)=pp(index(j),1:nx);
        elseif z(j)>z_alpha && z(j)<=z_beta
            z_beta=z(j);
            beta(1,1:nx)=pp(index(j),1:nx);
        elseif z(j)>z_beta && z(j)<=z_delta
            z_delta=z(j);
            delta(1,1:nx)=pp(index(j),1:nx);
        end
    end
    
    % Sorting the objective values and their corresponding solutions 
    for j=1:np-1
        for jj=j+1:np
            if z(jj)<z(j)
                c1=z(j);
                z(j)=z(jj);
                z(jj)=c1;
                c2=index(j);
                index(j)=index(jj);
                index(jj)=c2;
            end
        end
    end
   
    % Calculating the diversity index for each solution (Eq. (14))
    for j=1:np
        k=0;
        for jj=1:np
            if (jj~=j) 
                k=k+1;
                ss(j)=abs(z(j)-z(jj));
                if k==1
                    min_ss=ss(j);
                elseif ss(j)<min_ss
                    min_ss=ss(j);
                end
            end
        end
        if k==0
            div(j)=inf;
        else
            div(j)=min_ss;
        end
    end
    
    % Determining three most-diversified solutions in the current run
    div_alpha=-inf;
    div_beta=-inf;
    div_delta=-inf;
    for j=1:np
        if div(j)>=div_alpha
            div_alpha=div(j);
            alpha_div(1,1:nx)=pp(index(j),1:nx);
        elseif div(j)<div_alpha && div(j)>=div_beta
            div_beta=div(j);
            beta_div(1,1:nx)=pp(index(j),1:nx);
        elseif div(j)<div_beta && div(j)>=div_delta
            div_delta=div(j);
            delta_div(1,1:nx)=pp(index(j),1:nx);
        end
    end
   
    % Repositioning the first half of the best-fitted solutions around the three most diversified ones (Eqs. (15-17))
    for j=1:np
        if j<=nbest
            jjj=index(j);
            random=rand(1,1);
            if random<=(1/3)
                guide(1,1:nx)=alpha_div(1,1:nx)+sign(sign(rand(1,1)-0.5)+0.5)*((varmax(1,1:nx)-varmin(1,1:nx)).*rand(1,nx)+varmin(1,1:nx));
            elseif random>(1/3) && random<=(2/3)
                guide(1,1:nx)=beta_div(1,1:nx)+sign(sign(rand(1,1)-0.5)+0.5)*((varmax(1,1:nx)-varmin(1,1:nx)).*rand(1,nx)+varmin(1,1:nx));
            elseif random>(2/3) && random<=1
                guide(1,1:nx)=delta_div(1,1:nx)+sign(sign(rand(1,1)-0.5)+0.5)*((varmax(1,1:nx)-varmin(1,1:nx)).*rand(1,nx)+varmin(1,1:nx));
            end
            pp(jjj,1:nx)=guide;
        end
    end
    z_optimal(it)=z_alpha;
    optimal_pos(it,:)=alpha(1,:);

    % Saving the best-so-far objective value in the current run
    z_iter(it)=z_optimal(it);
end

% Save the final best solution and objective revealed upon the end of the optimization process
z_final=z_optimal(maxit);
pos_final(1:nx)=optimal_pos(maxit,1:nx);
end