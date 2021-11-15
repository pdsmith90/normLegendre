function [P]=normLegendre(order,t)
% recursion formula to evaluate 
% the fully normalized associated Legendre function
%
% https://link.springer.com/content/pdf/10.1007%2Fs00190-002-0216-2.pdf
%
% input: maximum order, l>=2
% input: t=vector
%
% output: matrix of coefficients, shifted by one, so P(1,1)=P_0,0
%         row-1 = l
%         column-1 = m
%         third dimension for each elemet of t


%how many input ts? (how long is 3rd dimension of coeff matrix)
k=length(t);

% initialize P matrix
P=zeros(order+1,order+1,k);

% calculate u=sin(t)
u=sin(t);

% initial values
% P_0,0
P(1,1,1:k)=ones(1,k);
% P_1,1
P(2,2,1:k)=sqrt(3).*u;

% !!!!
% P_1,0
%P(2,1,1:k)=sqrt(6).*cos(t);
P(2,1,1:k)=zeros(1,k);

% sectorals = diagonals, where m=l
for in=3:(order+1)
    % !!! actual n
    n=in-1;
    % here and below, have to use reshape because funky 3rd dimension
    % matrix element-wise multiplication
    % .*sqrt(1-t.^2) removed because not in paper but in lab prompt?
    P(in,in,1:k)=u.*sqrt((2.*n+1)./(2.*n)).*reshape(P(in-1,in-1,1:k),1,[]);
end


% iterate through l
% standard forward column recursion method
for in=3:(order+1)
    % !!! actual l
    l=in-1;
    % for each l, iterate through m for m<l
    for jn=1:(in-1)
        % !!! actual m
        m=jn-1;
        % dont have to save alpha and betas, 
        % just overwrite same variable for each iteration
        al=sqrt(((2.*l-1).*(2.*l+1))./((l-m).*(l+m)));
        be=sqrt(((2.*l+1).*(l+m-1).*(l-m-1))./((2.*l-3).*(l-m).*(l+m)));
        % compute P_l,m (remeber offset index)
        % again, funky reshapes
        P(in,jn,1:k)=al.*cos(t).*reshape(P(in-1,jn,1:k),1,[])-be.*reshape(P(in-2,jn,1:k),1,[]);   
    end
end

end
