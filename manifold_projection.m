function  y=manifold_projection(Wn)
[Nt K]=size(Wn);
%     y=sqrt(diag(diag(Wn*Wn')))^(-1)*Wn;
%         y=(sqrt(diag(diag(Wn'*Wn)))^(-1)*Wn')';
for jj=1:10,
    Wn=sqrt(diag(diag(Wn*Wn')))^(-1)*Wn;
    Wn=sqrt(Nt/K)*(sqrt(diag(diag(Wn'*Wn)))^(-1)*Wn')';
end
y=Wn;