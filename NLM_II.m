function [N_processed] = NLM_II(N,Ds,ds,h)
%% Non-local Means filter Accelerated by Integral Image:
%% Initialization:
% input the original data:
input=N;
% parameters:
[m,n]=size(input);
d=2*ds+1;% size of the comparing window.
D=2*Ds+1;% size of the searching window.
% h means the filter parameter, bigger h will filter more infprmation about
% texture, which can lead to the ambiguity of the porcessed result.
%%
% Enlarge the size of the input data N form n¡Ám to 
% (n+2Ds+2ds+1)¡Á(m+2Ds+2ds+1), and the new matrix called N_enlarged:
N_enlarged=padarray(input,[ds+Ds+1,ds+Ds+1],'symmetric','both');
%% Assign initial values to matrix that needed:
N_pocessing=zeros(m,n);
N_f=zeros(m,n);
maxweight=zeros(m,n);
Static_patch=N_enlarged(1+Ds:Ds+m+2*ds+1,1+Ds:Ds+n+2*ds+1);
%% Main loop:
for p=-Ds:Ds
    for q=-Ds:Ds       
        %Skip the current calculation when (p,q) equals (i,j):
        if(p==0&&q==0)
            continue;
        end
        % Calculate the integral image:
        Moving_patch=N_enlarged(1+Ds+p:Ds+m+2*ds+1+p,1+Ds+q:Ds+n+2*ds+1+q);
        difference=(Static_patch-Moving_patch).^2;
        Integral_image=cumsum(cumsum(difference,1),2);
        % Calculate the Euclidean distance:
        S=Integral_image(2*ds+2:m+2*ds+1,2*ds+2:n+2*ds+1)+Integral_image(1:m,1:n)-Integral_image(1:m,2*ds+2:n+2*ds+1)-Integral_image(2*ds+2:m+2*ds+1,1:n);
        S=S./(d^2); 
        % Calculate the weight:
        weight=exp(-S./(h*h));
        N_pocessing=N_pocessing+weight.*Moving_patch(ds+2:ds+1+m,ds+2:ds+1+n);
        N_f=N_f+weight;
        maxweight=max(maxweight,weight);
    end
end
N_pocessing=N_pocessing+maxweight.*Static_patch(ds+2:ds+1+m,ds+2:ds+1+n);
N_f=N_f+maxweight;
% Calculate the final result:
N_processed=N_pocessing./N_f;
end

