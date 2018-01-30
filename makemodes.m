% Chris
%matlabpool open local
%latsize = 8888;
%load h01.mat;
%I=speye(49152);
%k = 0.1570;
%M = I-(k*B);
%clear I k B
latsize

%%%%%%%%%%%%% 8888 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(latsize == 8888)
  % noise for Mx=z
 % z = z4noise(rank,1);
% load z2.mat; $ Chris
z = z2(:,1);
  % tolerance for eigenmdoes
  modetol = 1e-8;

%Chris -- new Arguments that make changes simplier
jenny = 200;
amy = 180; %changed from 160 to 180 TW 9/27/16

  % solve for modes in parallel
  spmd(2,2)
    if (labindex == 1)
      [x,d,vr]=gmresdrEIG(M,z,jenny,amy,1e-150,100);% Chris -- jenny/amy, changed parameters here as below, TW 9/28/17
    elseif (labindex == 2)
      [xh,gd,gvr,ghk,gvk]=hybridHF(M,z,jenny,amy,1e-150,100); % Chris,jenny/amy, changed rtol to 1e-150 from 1e-50 and cyclim 
    end                                                         % from 150 to 400 TW 9/28/17
  end

  % store mode information from various nodes
  clear x xh
  d   = d{1};   vr  = vr{1};
  gd  = gd{2};  gvr = gvr{2};
  ghk = ghk{2}; gvk = gvk{2};

  %%%%%%%%%%%%%% NON-HERMITIAN MODES %%%%%%%%%%%%%%%%
  % identify only acceptable Right modes
  res = zeros(amy,1); % Chris -- amy
  parfor ii=1:amy % Chris -- amy
    res(ii) = norm(M*vr(:,ii) - d(ii)*vr(:,ii));
  end

  % sort residuals
  [dsrt,indx] = sort(res);

  % find eligible modes
  ii = 0;
  while ( (ii < amy) && (dsrt(ii+1) < modetol) ) % Chris -- amy
    ii = ii + 1;
  end

  % remove inelegible modes
  d  = d(indx(1:ii));
  vr = vr(:,indx(1:ii));

  % sort modes by eig value magnitude
  [dsrt,indx] = sort(abs(d));
  d  = d(indx);
  vr = vr(:,indx);

  % find pairs by magnitude (should gaurentee even number of modes)
  clear indx;
  cnt = 0;
  for jj=1:ii-1
    for kk=(jj+1):ii
      if (abs(d(jj) - conj(d(kk))) < modetol)
        cnt = cnt + 1;
        indx(cnt) = jj;
        cnt = cnt + 1;
        indx(cnt) = kk;
% if we have 'if' command should we necessarily have 'else' command ?
      end
    end
  end
  d  = d(indx);
  vr = vr(:,indx); 
  ii = cnt;

  % make left eigen modes
  indx = 1:ii;
  for jj=1:2:ii
    tmp = indx(jj);
    indx(jj) = indx(jj+1);
    indx(jj+1) = tmp;
  end
  vl = gammaX(vr(:,indx),5);

%  resvl = zeros(1,size(vl,2));
 % parfor n = 1:size(vl,2)
  %  resvl(n) = norm(vl(:,n)'*M-d(n)*vl(:,n)');
 % end


  % normalize L/R modes, changed to 1/sqrt(tmp) by TW 9/18/17
 % for jj=1:ii
  %  tmp(jj) = vl(:,jj)' * vr(:,jj);
   % vr(:,jj) = (1/sqrt(tmp(jj))) * vr(:,jj);
   % vl(:,jj) = (1/sqrt(tmp(jj))) * vl(:,jj);
 % end

  % record number of acceptable modes
  neig = ii;

  %calculate eigenvalues of M*gamma5 using analytical soln
  %[heig, veig] = eigcalc(vl,vr,d,M);

  %%%%%%%%%%%%%%% HERMITIAN MODES %%%%%%%%%%%%
  % identify only acceptable Right modes
  res = zeros(amy,1); % Chris -- amy
  parfor ii=1:amy % Chris -- amy
    res(ii) = norm(M*gammaX(gvr(:,ii),5) - gd(ii)*gvr(:,ii));
  end

  %for ii = 1:amy
   % disp(res(ii))
  %end

  % sort residuals
  [dsrt,indx] = sort(res);

  % find eligible modes
  ii = 0;
  while ( (ii < amy) && (dsrt(ii+1) < 1e-6) ) % Chris -- amy, changed tolerance to accept more modes TW 9/25/17
    ii = ii + 1;
  end

  % remove inelegible modes
  gd  = gd(indx(1:ii));
  gvr = gvr(:,indx(1:ii));

  % sort by abs of eig value
  [dsrt,indx] = sort(abs(gd));
  gd  = gd(indx);
  gvr = gvr(:,indx);

  % record number of acceptable modes
  neigH = ii;

  % memory cleanup
  clear modetol z ii res dsrt indx cnt jj kk tmp
else
  return; 
end
% Chris -- stop clears, saves, and ending matlab
%clear M;
%save('eigen.mat')
%save('eigen.mat','vl','vr','gd','gvr')
% matlabpool close


























