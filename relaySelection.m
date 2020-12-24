function ind_relay=relaySelection(H,Pt)

[A_U,NT,N_user] = size(H) ;
N_legi=N_user-1;
h1=H(:,:,1);
h2=H(:,:,2);
P_common=Pt*0.5;
P_private_k=(Pt-P_common)/N_legi;
H_tot=[h1' h2'];%NT*2
[U2,~,~]=svd(H_tot);
hat_p_c=U2(:,1);

p_1=h1'/norm(h1)*sqrt(P_private_k);
p_2=h2'/norm(h2)*sqrt(P_private_k);
p_c=hat_p_c*sqrt(P_common); 

R_c1=log2(1+abs(h1*p_c)^2/(abs(h1*p_1)^2+abs(h1*p_2)^2+1));
R_c2=log2(1+abs(h2*p_c)^2/(abs(h2*p_1)^2+abs(h2*p_2)^2+1));

[~,ind_relay]=max([R_c1,R_c2]);

end