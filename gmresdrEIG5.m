function [x,th,evector, hk,vk,greal] = gmresdrEIG5(A,b,m,k,rtol,cyclim,fid) 
% gmresdr   solves systems of linear equations using a deflated restarted 
%           method. 
%
% Synopsis: 
%     [x,vk,hk,gdr] = gmresdr(A,b)
%     [x,vk,hk,gdr] = gmresdr(A,b,m,k)
%     [x,vk,hk,gdr] = gmresdr(A,b,m,k,rtol,cyclim)
%
% Input:
%   A      = an nxn non-hermitian matrix
%   b      = single right hand side
%   m      = size of Krylov Subspace (number of cycles before restarting)
%            Default = 10
%   k      = number of deflated eigenvalues/eigenvectors
%            Default = 5
%   rtol   = relative residual nececcary
%            Default = 1e-6
%   cyclim = maximum number of cycles to perform
%            Default = 100
%
% Output:
%   x   = solution vector to Ax=b
%   vk  = k approximate eigenvectors
%   hk  = k approcimate eigenvalues (located on diagonal)
%   gdr = residual at each iteration
%   
% Notes: In order to use gmresdr with gmresproj the relative residual
% should be run excessivly high to gain more accurate eigenvalue/vector
% information. This will lead to much faster convergence of gmresproj.
% TW changed all v's to vh for testing purposes 5/9/17

n = size(A,1);

if (nargin < 7)
  fid=1;
  if (nargin < 5)
    rtol = 1e-6;
    cyclim = 100;
    if (nargin < 3)
      m = 10;
      k = 5;
    end
  end
end

x = zeros(n,1);
cycle = 1;
j = 1;

rninit = norm(b(:,1));
rn = rninit;
r = b;
vn = norm(r);
vh(:,1) = r/vn;
c(1,1) = vn;
c(2:m+1,1) = zeros(m,1);

while ( (cycle <= cyclim) )
%while ( (rn/rninit > rtol) && (cycle <= cyclim) )
  while ( j <= m ) 
    wv =A* gamma5(vh(:,j),n,1);
    vnf = norm(wv);

    for i = 1:j
      h(i,j) = vh(:,i)'*wv;
      wv = wv - h(i,j) * vh(:,i);
    end
    vn = norm(wv);

    %--------reorthogonalization section-------%
    if( vn < 1.1*vnf )
      for i = 1:j
        dot = vh(:,i)'*wv;
        wv = wv - dot * vh(:,i);
        h(i,j) = h(i,j) + dot;
      end
      vn = norm(wv);
    end
    %--------------------------------------------------%

    h(j+1,j) = vn;
    vh(:,j+1) = wv/h(j+1,j);

  % output data per cycle
  % disp(['iteration: ',num2str(j)]);

    j = j + 1;
  end
  
  j = 1;
  
 % for i = 1:m
  %   d(i,1) = v(:,i)'*(gamma5(xvec(j,:)',n,1)-x(:,1));%added to back compute the vector d
 % end

  d(1:m,1) = h(1:m+1,1:m) \ c(1:m+1) ;
  srv(1:m+1,1) = c(1:m+1)-(h(1:m+1,1:m)*d(1:m,1));

  %----Set up and solve linear equations.-----%
  x(:,1) = x(:,1) + vh(:,1:m)*d(1:m);
 % x(:,1) = xvec(j,:)';
  r = vh(:,1:m+1)*srv;
  rn = norm(r);

  j = j + 1;

  hh = h(1:m,1:m);
  em = zeros(m,1);
  em(m) = 1;
  ff = hh' \ em;

  hh(:,m) = hh(:,m) + h(m+1,m)^2 * ff;
  hh=full(hh); [g,dd] = eig(hh,'nobalance'); 
% no balance needed!!!-WW; 03/25/2014!

  dd = diag(dd);
  dabs = abs(dd);
  [dabs,ind] = sort(dabs);
  th = dd(ind);
  gg = g(:,ind(1:k));
 
  for i=1:k
    rho(i) = gg(:,i)'*h(1:m,:)*gg(:,i);
    tv = h(1:m,:)*gg(:,i)-rho(i)*gg(:,i);
    tvn = norm(tv);
    rna(cycle,i) = sqrt( tvn*tvn+ abs(h(m+1,m))^2*abs(gg(m,i))^2 );
    tha(cycle,i) = th(i);
    rhoa(cycle,i) = rho(i);
  end

  [rnasort,ind] = sort(rna(cycle,:));
  rna(cycle,1:k)=rna(cycle,ind(1:k));
  gg = gg(:,ind(1:k));
  th = th(ind);

  greal = gg;

  greal(m+1,1:k) = zeros(1,k);
  greal(:,k+1) = srv;

  [gon,rr] = qr(greal(:,1:k+1),0);
  hnew = gon'*h*gon(1:m,1:k);
  h(k+1,:) = zeros(1,m);

  j = 1;
  rtolev = 1e-11;
  rtolev = 1e-14;
  while (j <= k && rna(cycle,j) <= rtolev)
    hnew(j+1:k+1,j) = zeros(k-j+1,1);
    j = j + 1;
  end



  evector = vh(:,1:m) * greal(1:m,1:k);


  h(1:k+1,1:k) = hnew;

  c(1:k+1,1) = gon(:,1:k+1)'*srv(:,1);
  c(k+2:m+1,1) = zeros(m-k,1);

  work = vh*gon;
  vh(:,1:k+1) = work;

  %section for just reorthog. one vector, v_{k+1}
  for i = 1:k
    dot = vh(:,i)'*vh(:,k+1) ;
    vh(:,k+1) = vh(:,k+1) - dot * vh(:,i);
  end

  vh(:,k+1) = vh(:,k+1)/norm(vh(:,k+1));

  % output data per cycle
  %fprintf(fid,'Cycle: %d  Rel Res Norm: %12.8g\n',cycle,rn/rninit);

  j = k + 1;
  cycle = cycle + 1;
end

hk = h(1:k+1,1:k);
vk=vh(:,1:k+1);

%for i=1:k
% dispeig = ['value of eig ', i,' ', num2str(th(i))];
 %disp(dispeig);
%end

%semilogy(rna)                                                               

if (rn/rninit < rtol) && (cycle-1 <= cyclim)
  fprintf(fid,'gmresdrEIG5(%d,%d) converged in %d cycles with a relative residual of %12.8g\n', ...
          m,k,cycle-1,rn/rninit);
else
  fprintf(fid,'gmresdrEIG5(%d,%d) stoped after %d cycles without converging to the desired tolerence of %12.8g\n', ...
          m,k,cycle-1,rn/rninit);
  fprintf(fid,'a relative residual of %12.8g has been reached\n',rn/rninit);
end

return
