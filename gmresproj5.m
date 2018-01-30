function [x,gdr] = gmresproj5(A,b,m,k, vk, hk,rtol,cyclim,fid)
% gmresproj solves systems of linear equations using a deflated restarted 
%           method, where the approximate eigenvector/values have been
%           calculated already by gmresdr()
%
% Synopsis: 
%     [x,gdr] = gmresproj(A,b,m,k,vk,hk)
%     [x,gdr] = gmresproj(A,b,m,k,vk,hk,rtol,cyclim)
%
% Input:
%   A      = an nxn non-hermitian matrix
%   b      = single right hand side
%   m      = size of Krylov Subspace (number of cycles before restarting)
%   k      = number of deflated eigenvalues/eigenvectors
%   vk     = k approximate eigenvectors (from gmresdr)
%   hk     = k approximate eigenvalues (located on diagonal) (from gmresdr)
%   rtol   = relative residual nececcary
%            Default = 1e-6
%   cyclim = maximum number of cycles to perform
%            Default = 100
%
% Output:
%   x   = solution vector to Ax=b
%   gdr = residual at each iteration
%   
% Notes: In order to use gmresproj, gmresdr must have already been run. For
% better convergence gmresdr should be run with small residual to gather
% very accurate approximate eigenvector/values. 

if (nargin < 9)
  fid = 1;
  if (nargin < 7)
    rtol = 1e-6;
    cyclim = 100;
  end
end


mvp=0;

x=zeros(size(A,1),1);
h(1:m+1,1:m) = zeros(m+1,m);
%for irhs = 2:nrhs
cycle = 1;
j = 1;
cmult = 1;
rninit = norm(b(:));
rn = rninit;
r = b; %minv(n,b(:));

while ( (rn/rninit > rtol) & (cycle <= cyclim) ) 
%
%  Projection section (keep res's parallel for diff. shifts)
%
    iproj = 1;
	ifreq = 1;
    if( iproj ~= 0 ) 
    if( cycle ==  floor(cycle/ifreq) * ifreq )
    c(1:k+1,1) = vk'*r;
	eyebar = eye(k+1,k);
    eyenotbar = eye(k,k);
    %
    %  Question on which projection should be used here.  Must
    %  decide!!!!!
  	d(1:k,1) = (hk) \ c(1:k+1,1) ;
    srv(1:k+1,1) = c(1:k+1,1)-(hk)*d(1:k,1);
	srvhere = srv(1:k+1,1);

		x(:) = x(:) + vk(:,1:k)*d(1:k);
		wv = b(:) -A* gamma5(x,size(A,1),1); %op(n,x(:,irhs)); %-may remove this and next line%
% shifted A out of gamma5		
	r = wv;
end  %if involving ifreq
end  %if involving iproj
%
%  end projection section
%
vn = norm(r);
v(:,1) = r/vn;
  
c(1,1) = vn;
c(2:m+1,1) = zeros(m,1);

while ( (j <= m) ) 
    wv =A* gamma5(v(:,j),size(A,1),1); %op(n,v(:,j));
%BS    wv = gamma5(A*v(:,j),size(A,1),1); %op(n,v(:,j));
%BS shifted A out of gamma5   
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
         % disp( 'did a reorthog')
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
    eyenotbar = eye(j,j);
%subsection for res.norm at every iteration
 	d(1:j,1) = (h(1:j+1,1:j)) \ c(1:j+1) ;
    srv(1:j+1,1) = c(1:j+1)-(h(1:j+1,1:j))*d(1:j,1);

    gdr(mvp,1) = norm(srv(1:j+1,1));
    j = j + 1;
end  %while

  %----Set up and solve linear equations.-----%
%  Should already be done:  d = h \ c;

	x(:) = x(:) + v(:,1:m)*d(1:m);
	wv = b(:) -A* gamma5(x,size(A,1),1); %op(n,x(:,irhs)); %-may remove this and next line%
% BS shift A outside of gamma5        wv = b(:) - gamma5(A*x,size(A,1),1); %op(n,x(:,irhs)); %-may remove this
	% BS shifted A out of gamma5 rnale(cycle) = norm(wv);
% Note could change here to compute the true residual:
r = v(:,1:m+1)*srv(1:m+1,1);
rn = norm(r);


%disp(['Cycle: ',num2str(cycle), ' Rel Res Norm: ', num2str(rn/rninit)])


v(:,1) = r/rn;
c(1,1) = rn;
c(2:m+1,1) = zeros(m,1);

j = 1;
cycle = cycle + 1;
end  %while


if (rn/rninit < rtol) && (cycle-1 <= cyclim)
  fprintf(fid,'gmresproj5(%d,%d) converged in %d cycles with a relative residual of %12.8g\n', ...
          m,k,cycle-1,rn/rninit);
else
  fprintf(fid,'gmresproj5(%d,%d) stoped after %d cycles without converging to the desired tolerence of %12.8g\n', ...
          m,k,cycle-1,rn/rninit);
  fprintf(fid,'a relative residual of %12.8g has been reached\n',rn/rninit);
end


return 
