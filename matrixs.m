function [varargout]= matrixs(p)

Ap=zeros(p.Np-1);
An=zeros(p.Nn-1);
  
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


Bp= zeros(p.Np-1,2);
Bp(end,end)=(p.Np)/(p.Np-1) * sp;

Bn= zeros(p.Nn-1,2);
Bn(end,end)=(p.Nn)/(p.Nn-1) * sn;

N1 = zeros(2,p.Nn-1);
N1(1,1) = 4;
N1(1,2) = -1;
N1(2,end) = -4;
N1(2,end-1) = 1;

N2 = diag([-3,3]);

N3_n = [0; -(2*p.delta_n)/(p.Ds_n)];
N3_p = [0; -(2*p.delta_p)/(p.Ds_p)];

% N3_n = [0; (p.delta_n)/(p.Ds_n) ];
% N3_p = [0; -(p.delta_p)/(p.Ds_p)];



%A,B matrices for each electrode
A_n = An - Bn*(N2\N1);

A_p = Ap - Bp*(N2\N1);



B_n = Bn*(N2\N3_n);
B_p = Bp*(N2\N3_p);

%C,D matrices for each electrode
C_n = -[0,1]*(N2\N1);
C_p = -[0,1]*(N2\N1);

D_n = [0,1]*(N2\N3_n);
D_p = [0,1]*(N2\N3_p);


varargout{1}=A_n;
varargout{2}=A_p;
varargout{3}=B_n;
varargout{4}=B_p;
varargout{5}=C_n;
varargout{6}=C_p;
varargout{7}=D_n;
varargout{8}=D_p;


end
