function [varargout]=matrixe(p)

for i= 1:p.Nxn+p.Nxs
     
    A(i,i)=-2*(p.D_en_eff/(p.epsilon_e_n*p.Del_xn^2));
    for k= p.Nxn:p.Nxn-1+p.Nxs-1  %Nxn:Nxs+Nxn
        A(k,k)=-2*(p.D_es_eff/(p.epsilon_e_s*p.Del_xs^2));
    end
    for k=p.Nxn-1+p.Nxs:p.Nx-3  %Nxs+Nxn:Nx
        A(k,k)=-2*(p.D_ep_eff/(p.epsilon_e_p*p.Del_xp^2));
    end

    if (i~=1)
    A(i,i-1)=(p.D_en_eff/(p.epsilon_e_n*p.Del_xn^2));
    for k=p.Nxn:p.Nxn-1+p.Nxs-1  %Nxn:Nxs+Nxn
        A(k,k-1)=(p.D_es_eff/(p.epsilon_e_s*p.Del_xs^2));
    end
    for k=p.Nxn-1+p.Nxs:p.Nx-3  %Nxs+Nxn:Nx
        A(k,k-1)=(p.D_ep_eff/(p.epsilon_e_p*p.Del_xp^2));
    end
    
    end
    A(i,i+1)=(p.D_en_eff/(p.epsilon_e_n*p.Del_xn^2));
    for k=p.Nxn:p.Nxn-1+p.Nxs-1  %Nxn:Nxs+Nxn
        A(k,k+1)=(p.D_es_eff/(p.epsilon_e_s*p.Del_xs^2));
    end
     for k=p.Nxn-1+p.Nxs:p.Nx-4  %Nxs+Nxn:Nx-1
        A(k,k+1)=(p.D_ep_eff/(p.epsilon_e_p*p.Del_xp^2));
     end
end    

A(1,1)=-p.D_en_eff/(p.epsilon_e_n*p.Del_xn^2); 
A(p.Nx-3,p.Nx-3)=-p.D_ep_eff/(p.epsilon_e_p*p.Del_xp^2);


varargout{1}=A;



end
