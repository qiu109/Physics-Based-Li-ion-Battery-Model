function [varargout]= matrixs(p)


Ap=zeros(p.Np-1);
An=zeros(p.Nn-1);
Bp= zeros(p.Np-1,1);  
Bn= zeros(p.Nn-1,1);
%% Positive and Negative solid particle matrices

 sp=p.Ds_p/(p.delta_p^2);
 sn=p.Ds_n/(p.delta_n^2);
 
for k=1:p.Np-1
    
    %diagonal
    Ap(k,k)=-2*sp; 
    An(k,k)=-2*sn; 
    
    %superdiagonal
    if (k~= (p.Np-1))
    Ap(k,k+1)= (k+1)/k.*sp;  
    An(k,k+1)= (k+1)/k.*sn;
    end
    
    %subdiagonal
    if (k~= 1)
    Ap(k,k-1)= (k-1)/k.*sp;  
    An(k,k-1)= (k-1)/k.*sn;
    end
end


%1st order boundary discretization

% An(p.Nn-1,p.Nn-1)=(p.Ds_n/((p.Nn-1)*p.delta_n*(p.delta_n^2)))*(p.delta_n-(p.Nn-1)*p.delta_n); %equivalently sn*(2-p.Nn)/(p.Nn-1);
% Ap(p.Np-1,p.Np-1)=(p.Ds_p/((p.Np-1)*p.delta_p*(p.delta_p^2)))*(p.delta_p-(p.Np-1)*p.delta_p); %equivalently sp*(2-p.Np)/(p.Np-1);
% 
% Bp(end,end)=  -(p.Np)/(p.Np-1)   /p.delta_p;
% Bn(end,end)=  -(p.Nn)/(p.Nn-1)   /p.delta_n;


%2nd order boundary discretization

An(p.Nn-1,p.Nn-1)=(p.Ds_n*(6-2*p.Nn))/(3*(p.Nn-1)*p.delta_n^2);
An(p.Nn-1,p.Nn-2)=-(p.Ds_n*(6-2*p.Nn))/(3*(p.Nn-1)*p.delta_n^2);

Ap(p.Np-1,p.Np-1)=(p.Ds_p*(6-2*p.Np))/(3*(p.Np-1)*p.delta_p^2);
Ap(p.Np-1,p.Np-2)=-(p.Ds_p*(6-2*p.Np))/(3*(p.Np-1)*p.delta_p^2);

Bp(end,end)=  -(p.Np)/(p.Np-1) *  (2)/(3*p.delta_p);
Bn(end,end)=  -(p.Nn)/(p.Nn-1) *  (2)/(3*p.delta_n);

varargout{1}=Ap;
varargout{2}=An;
varargout{3}=Bn;
varargout{4}=Bp;

end
