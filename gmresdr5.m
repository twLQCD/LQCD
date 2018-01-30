function [x,vk,hk,gdr] = gmresdr5(A,b,m,k,rtol,cyclim,fid)
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

n = size(A,1);

if (nargin < 7)
  fid = 1;
  if (nargin < 5)
    rtol = 1e-6;
    cyclim = 100;
    if (nargin < 3)
      m = 10;
      k = 5;
    end
  end
end


kperm = k;
x = zeros(n,1);
cycle = 1;
mvp = 0;
j = 1;  % is this j the culprit? ask carl BS
rninit = norm(b(:,1));
rn = rninit;
r = b;
vn = norm(r);
v(:,1) = r/vn;
c(1,1) = vn;
c(2:m+1,1) = zeros(m,1);
cmult(1) = ones(1,1);
while ( (rn/rninit > rtol) && (cycle <= cyclim) ) 
while ( (j <= m) )
    wv = A*gamma5(v(:,j),n,1); % op(n,v(:,j)) ;% shifted A out of gamma5
    f = wv; %minv(n,wv); 
    mvp = mvp + 1;
    vnf = norm(f);

    for i = 1:j
      h(i,j) = v(:,i)'*f;
      f = f - h(i,j) * v(:,i);
  end  %for
    vn = norm(f);
%--------------------------------------------------%
%--------reorthogonalization section-------%
    if( vn < 0*vnf )
          fprintf(fid, 'did a reorthog')
      for i = 1:j
        dot = v(:,i)'*f;
        f = f - dot * v(:,i);
        h(i,j) = h(i,j) + dot;
      end  %for
      vn = norm(f);
    end  %if
%--------------------------------------------------%
    h(j+1,j) = vn;
    v(:,j+1) = f/h(j+1,j);
	eyebar = eye(j+1,j);
%subsection for res.norm at every iteration
	d(1:j,1) = (h(1:j+1,1:j)) \ c(1:j+1) ;
    srv(1:j+1,1) = c(1:j+1)-(h(1:j+1,1:j))*d(1:j,1);
    gdr(mvp,1) = norm(srv(1:j+1,1));


    j = j + 1;
end  %while

  %----Set up and solve linear equations.-----%
%  Should already be done:  d = h \ c;

	x(:,1) = x(:,1) + v(:,1:m)*d(1:m);
	wv = b(:,1) -A* gamma5(x,n,1); % op(n,x(:,1)); %-may remove this and next line% BS shifted A out of gamma5
	rnale(cycle) = norm(wv);
r = v(:,1:m+1)*srv;
rn = norm(r);


%disp(['Cycle: ',num2str(cycle), ' Rel Res Norm: ', num2str(rn/rninit)])

%
  hh = h(1:m,1:m);
  em = zeros(m,1);
  em(m) = 1;
  ff = hh' \ em;
%
%FOM section
%Comment next line for FOM instead of GMRES:
  hh(:,m) = hh(:,m) + h(m+1,m)^2 * ff;
  hh=full(hh); [g,dd] = eig(hh,'nobalance'); 
% no balance needed!!!-WW; 03/25/2014!
  ddd = diag(dd);
  dabs = abs(ddd);
  [dabss,ind] = sort(dabs);
  th = ddd(ind);
  gg = g(:,ind(1:k));
  for i=1:k
    rho(i) = gg(:,i)'*h(1:m,:)*gg(:,i);
    tv = h(1:m,:)*gg(:,i)-rho(i)*gg(:,i);
    tvn = norm(tv);
    rna(cycle,i) = sqrt( tvn*tvn+ abs(h(m+1,m))^2*abs(gg(m,i))^2 );
    tha(cycle,i) = th(i);
    rhoa(cycle,i) = rho(i);
  end  %for

    greal = gg;

  greal(m+1,1:k) = zeros(1,k);
  greal(:,k+1) = srv;
  [gon,rr] = qr(greal(:,1:k+1),0);
  hnew = gon'*h*gon(1:m,1:k);
  h(kperm+1,:) = zeros(1,m);
  h(1:k+1,1:k) = hnew;
  c(1:k+1,1) = gon(:,1:k+1)'*srv(:,1);
  c(k+2:m+1,1) = zeros(m-k,1);

  work = v*gon;
  v(:,1:k+1) = work;
  %section for just reorthog. one vector, v_{k+1}
  for i = 1:k
    dot = v(:,i)'*v(:,k+1) ;
    v(:,k+1) = v(:,k+1) - dot * v(:,i);
  end  %for
  v(:,k+1) = v(:,k+1)/norm(v(:,k+1));
  j = k + 1;
  kold = k;
  k = kperm;

  % output data per cycle
  %disp(['Cycle: ',num2str(cycle),' Rel Res Norm: ',num2str(rn/rninit)]);
  %fprintf(fid,'Cycle: %d  Rel Res Norm: %12.8g\n',cycle,rn/rninit);

  cycle = cycle + 1;
end  %while

k = kold;
vk = v(:,1:k+1);
hk = h(1:k+1,1:k);

if (rn/rninit < rtol) && (cycle-1 <= cyclim)
  fprintf(fid,'gmresdr5(%d,%d) converged in %d cycles with a relative residual of %12.8g\n', ...
          m,k,cycle-1,rn/rninit);
else
  fprintf(fid,'gmresdr5(%d,%d) stoped after %d cycles without converging to the desired tolerence of %12.8g\n', ...
          m,k,cycle-1,rn/rninit);
  fprintf(fid,'a relative residual of %12.8g has been reached\n',rn/rninit);
end
return
