
% --- clear data holders ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear M gM z gx x gd gvr d vr vl invxi xh 
clear nsg_raw esg_raw ps_raw ns_raw es_raw esps_raw esgps_raw
clear nsAVG nsgAVG esAVG esgAVG psAVG espsAVG esgpsAVG
clear nsERR nsgERR esAVG esgERR psERR espsERR esgpsERR
clear neig neigH

% --- prep initial info ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kB = k*B;
I = speye(rank);
M = I - kB;
M = sparse(M);
clear B
latsize
% --- Holder for noise and solution vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntrials = 1;
% nrhs = 200; % Chris -- defined in evscript
 % BS 5/23/014 znoise = zeros(rank,ntrials * nrhs);
 % BS 5/23/014 gxsol  = zeros(rank,ntrials * nrhs);
 % BS 5/23/014 xsol   = zeros(rank,ntrials * nrhs);

% --- EIG INFO --%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % BS 5/23/014 makemodes;
makemodes; % Chris
% load eigen.mat; % Chris
gdinv = inv(diag(gd));
dinv = inv(diag(d));

% --- GMRES INFO ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rtolDR=1e-20;
rtolPj=1e-14;
cyclemax=300;
mDR=80; kDR=40;

% --- Solve 1st rhs (results are thrown away) ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = z2(:,1);
tic;
[gx,gvk,gdk] = gmresdr5(M,z,mDR,kDR,rtolDR,cyclemax);
toc;

% Generate 1/xi for ESPS method

invxi = zeros(neig,neig,nsub+1);
Mp = vr;
invxi(:,:,1) = diag(diag(vl' * Mp));
for ksub = 2:nsub+1
  Mp = vr + (kB * Mp);
  invxi(:,:,ksub) = diag(diag(vl' * Mp));
end

% Generate 1/gxi for gESPS method
invgxi = zeros(neigH,neigH,nsub+1);
ggvr = gvr;
Mp = ggvr;
% BS 11/30/2014 this line replaced by line below it invgxi(:,:,1) = diag(diag(gvr' * Mp));
invgxi(:,:,1) = diag(diag(gvr' * gammaX(Mp,5)));

for ksub = 2:nsub+1
  Mp = ggvr + (kB * Mp);
  %  BS 11/30/2014 this line replaced by line below it invgxi(:,:,ksub) = diag(diag(gvr' * Mp));
invgxi(:,:,ksub) = diag(diag(gvr' * gammaX(Mp,5)));
end
clear Mp ggvr
% invgxi % Chris 

%For invgxipoly
% load z2.mat; 
% load h01.mat;
%I=speye(49152);
%k = 0.1570;
%M = I-(k*B);
v = zeros(49152,1);
v(:,1)= z2(:,1);

p = 8; % BS i think p is equivalent to ksub
n = 49152;
for j=1:p
       v(:,j+1) = M * v(:,j);
   end
lsmat = v(:,2:p+1)'*v(:,2:p+1);
cls = v(:,2:p+1)'*v(:,1);
  dls = lsmat\cls;
  A = dls





invgxipoly = zeros(neigH,neigH,nsub+1);
ggvr = gvr;
Mp = A(1,1)*ggvr;
temp=gvr;




invgxipoly(:,:,1) = diag(diag(gvr' * gammaX(Mp,5)));


for ksub = 2:nsub+1
temp=M*temp;
Mp=Mp+A(ksub,1)*temp;
invgxipoly(:,:,ksub) = diag(diag(gvr' * gammaX(Mp,5)));

end





% --- Loop for various trials ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iruns=1:ntrials
  disp([' ---------- Run: ',num2str(iruns)]);

  % --- noise ---

%%%  z = z4noise(rank,nrhs);
%  load z2.mat; % Chris -- we now generate z2 in evscript
  znoise(:,(iruns-1)*nrhs+1:iruns*nrhs) = z2;
  z = z2;
  % --- Prep ES data with noise
  invgd_gvrdag_z = gdinv * (gvr' * z);
  invd_vldag_z = dinv * (vl' * z);

  % --- Prep PS data with noise
  Minvpert_z = zeros(rank,nrhs,nsub+1);
  Minvpert_z(:,:,1) = z;
  for ksub=1:nsub
    Minvpert_z(:,:,ksub+1) = z + (kB * Minvpert_z(:,:,ksub));
  end


% Now comes the part for polynomial method BS  4/3/2015
% Value of A was initially calculated elsewhere
%A = [   6.5614 - 0.0151i
 %     -21.0752 + 0.0696i
  %     40.3994 - 0.1562i
   %   -48.9310 + 0.2083i
    %   37.7279 - 0.1729i
     % -17.9445 + 0.0880i
      %  4.7977 - 0.0252i
       %-0.5516 + 0.0031i];
% Determine coefficients
% Chris -- comments out section mini-version
%load h01.mat
%I=speye(49152);
%k = 0.1570;
%M = I-(k*B);
% end Chris comment storm

% Chris -- changed from 
% v = zeros(rank,1);
% v(:,1)=z2(:,1);
A = zeros(rank,1);
A = z2(:,1);
v(:,1)= A;

p = 8; % BS i think p is equivalent to ksub
n = rank;
for j=1:p
       v(:,j+1) = M * v(:,j);
   end
lsmat = v(:,2:p+1)'*v(:,2:p+1);  
cls = v(:,2:p+1)'*v(:,1);
  dls = lsmat\cls;
  A = dls

 Minvper_z = zeros(rank,nrhs,nsub+1);
   Minvper_z(:,:,1) = A(1,:)*z;
tryn = zeros(rank,nrhs,nsub+1);
tryn(:,:,1) = z;

  for ksub=1:nsub

    tryn(:,:,ksub+1) = (M * tryn(:,:,ksub));

    Minvper_z(:,:,ksub+1) =  Minvper_z(:,:,ksub)+A(ksub+1,:)* tryn(:,:,ksub+1);
  end




 % --- Prep ESPS  and gESPS data with noise
  invxi_vldag_z = zeros(neig,nrhs,nsub+1);
  invgxi_gvrdag_z = zeros(neigH,nrhs,nsub+1);
  invgxipoly_gvrdag_z = zeros(neigH,nrhs,nsub+1);
  for ksub=1:nsub+1
    invxi_vldag_z(:,:,ksub) = invxi(:,:,ksub) * (vl' * z);
    invgxi_gvrdag_z(:,:,ksub) = invgxi(:,:,ksub) * (gvr' * z);
    invgxipoly_gvrdag_z(:,:,ksub) = invgxipoly(:,:,ksub) * (gvr' * z);
  end

  % --- Solution Vectors ---
  gx = zeros(rank,nrhs);
  parfor rhs=1:nrhs
    tic;
    disp(['rhs: ',num2str(rhs)]);
    gx(:,rhs) = gmresproj5(M,z(:,rhs),mDR,kDR,gvk,gdk,rtolPj,cyclemax);
    toc;
  end
  gxsol(:,(iruns-1)*nrhs+1:iruns*nrhs) = gx;


  % --- Run RHS ---

  nsg_raw   = zeros(nrhs,5);
  es_raw    = zeros(nrhs,neig,5);
  esg_raw   = zeros(nrhs,neigH,5);
  ps_raw    = zeros(nrhs,nsub+1,5);
  poly_raw  = zeros(nrhs,nsub+1,5);
  esps_raw  = zeros(nrhs,neig,nsub+1,5);
  esgps_raw = zeros(nrhs,neigH,nsub+1,5);
  JaveHFES  = zeros(nrhs,neigH,5);
  JaveHFESPS = zeros(nrhs,neigH,nsub+1,5);
  esgpoly_raw= zeros(nrhs,neigH,nsub+1,5);
  esgpolylatest_raw= zeros(nrhs,neigH,nsub+1,5);
  for rhs=1:nrhs
    tic; 
    disp([' ----- Rhs: ',num2str(rhs)]);

    % --- NS ---
    for iop=1:5
      if (iop == 5)  % scalar
        nsg_raw(rhs,iop) = z(:,rhs)'*gammaX(gx(:,rhs),5);
      else  % gamma_mu
        nsg_raw(rhs,iop) = z(:,rhs)'*gammaX(gammaX(gx(:,rhs),5),iop);
 % BD 11/14/2014 interchanged position of iop and 5
      end
    end

    % --- ESg (gamma5) ---
    xp = zeros(rank,1);
    for ieig=1:neigH
      xp = xp + invgd_gvrdag_z(ieig,rhs)*gvr(:,ieig);
      for iop = 1:5
        if (iop == 5) % scalar
          esg_raw(rhs,ieig,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(xp,5);
          JaveHFES(rhs,ieig,iop) = z(:,rhs)'*gammaX(xp,5);
        else % gamma_mu
          esg_raw(rhs,ieig,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(gammaX(xp,5),iop);
          JaveHFES(rhs,ieig,iop) =  z(:,rhs)'*gammaX(gammaX(xp,5),iop);
% esg_raw(rhs,ieig,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(gammaX(xp,iop),5);
 % BD 11/14/2014 interchanged position of iop and 5
       
       end
      end
    end

    % --- ES ---
    xp = zeros(rank,1);
    for ieig=1:neig
      xp = xp + invd_vldag_z(ieig,rhs)*vr(:,ieig);
      for iop = 1:5
        if (iop == 5) % scalar
          es_raw(rhs,ieig,iop) =  nsg_raw(rhs,iop) - z(:,rhs)'*xp;
        else % gamma_mu
          es_raw(rhs,ieig,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(xp,iop);
        end
      end
    end

    % --- PS ---
    for ksub=1:nsub+1
      for iop=1:5
        if (iop == 5) % scalar
          ps_raw(rhs,ksub,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*Minvpert_z(:,rhs,ksub); 
% is something missing in it?      
  else  % gamma_mu
          ps_raw(rhs,ksub,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(Minvpert_z(:,rhs,ksub),iop);
% --- clear data holders ---
% --- clear data holders ---
% --- clear data holders ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
      end
    end


% Now is the polynomial part BS 4/3/2015
  for ksub=1:nsub+1
      for iop=1:5
        if (iop == 5) % scalar
          poly_raw(rhs,ksub,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*Minvper_z(:,rhs,ksub);
% is something missing in it?
  else  % gamma_mu
          poly_raw(rhs,ksub,iop) = nsg_raw(rhs,iop) - z(:,rhs)'*gammaX(Minvper_z(:,rhs,ksub),iop);
% --- clear data holders ---
% --- clear data holders ---
% --- clear data holders ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
      end
    end









    % --- ESPS ---
    xp = zeros(rank,nsub+1);
    for ieig = 1:neig
      for ksub = 1:nsub+1
        xp(:,ksub) = xp(:,ksub) + invxi_vldag_z(ieig,rhs,ksub)*vr(:,ieig);
        for iop=1:5
          if (iop == 5) % scalar
            esps_raw(rhs,ieig,ksub,iop) = es_raw(rhs,ieig,iop) + ps_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * xp(:,ksub));
          else
            esps_raw(rhs,ieig,ksub,iop) = es_raw(rhs,ieig,iop) + ps_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(xp(:,ksub),iop));
         % should ksub and iop be interchanged?
           end
        end
      end
    end

    % --- ESgPS ---
    xp = zeros(rank,nsub+1);
    for ieig = 1:neigH
      for ksub = 1:nsub+1
        xp(:,ksub) = xp(:,ksub) + invgxi_gvrdag_z(ieig,rhs,ksub)*gvr(:,ieig);
        for iop=1:5
          if (iop == 5) % scalar
            esgps_raw(rhs,ieig,ksub,iop) = esg_raw(rhs,ieig,iop) + ps_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(xp(:,ksub),5));
            JaveHFESPS(rhs,ieig,ksub,iop) = (z(:,rhs)' * gammaX(xp(:,ksub),5));
          else
            esgps_raw(rhs,ieig,ksub,iop) = esg_raw(rhs,ieig,iop) + ps_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(gammaX(xp(:,ksub),5),iop));
            JaveHFESPS(rhs,ieig,ksub,iop) =  (z(:,rhs)' * gammaX(gammaX(xp(:,ksub),5),iop));
          end
        end
      end
    end


% This is the part where i worked on HFPOLY BS 5/6/2015

    xp = zeros(rank,nsub+1);
    for ieig = 1:neigH
      for ksub = 1:nsub+1
        xp(:,ksub) = xp(:,ksub) + invgxi_gvrdag_z(ieig,rhs,ksub)*gvr(:,ieig);
        for iop=1:5
          if (iop == 5) % scalar
            esgpoly_raw(rhs,ieig,ksub,iop) = esg_raw(rhs,ieig,iop) + poly_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(xp(:,ksub),5));
          else
            esgpoly_raw(rhs,ieig,ksub,iop) = esg_raw(rhs,ieig,iop) + poly_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(gammaX(xp(:,ksub),5),iop));
          end
        end
      end
    end




% This is the part where i worked on HFPOLYLATEST BS 5/6/2015

    xp = zeros(rank,nsub+1);
    for ieig = 1:neigH
      for ksub = 1:nsub+1
        xp(:,ksub) = xp(:,ksub) + invgxipoly_gvrdag_z(ieig,rhs,ksub)*gvr(:,ieig);
        for iop=1:5
          if (iop == 5) % scalar
            esgpolylatest_raw(rhs,ieig,ksub,iop) = esg_raw(rhs,ieig,iop) + poly_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(xp(:,ksub),5));
          else
            esgpolylatest_raw(rhs,ieig,ksub,iop) = esg_raw(rhs,ieig,iop) + poly_raw(rhs,ksub,iop) - nsg_raw(rhs,iop) + (z(:,rhs)' * gammaX(gammaX(xp(:,ksub),5),iop));
          end
        end
      end
    end














    toc;
  end


 %start correcting for the real and imaginary results BS 12/8/2015
            for iop = 1:5
                if (iop ==5)
                   esgpolylatest_raw(:,:,:,iop) = real(esgpolylatest_raw(:,:,:,iop));
                   esgpoly_raw(:,:,:,iop) = real(esgpoly_raw(:,:,:,iop));
                   esgps_raw(:,:,:,iop)   = real(esgps_raw(:,:,:,iop));
                   esps_raw(:,:,:,iop)    = real(esps_raw(:,:,:,iop));
                   poly_raw(:,:,iop)      = real(poly_raw(:,:,iop));
                   ps_raw(:,:,iop)        = real(ps_raw(:,:,iop));
                   es_raw(:,:,iop)        = real(es_raw(:,:,iop));
                   esg_raw(:,:,iop)       = real(esg_raw(:,:,iop));
                   nsg_raw(:,iop)         = real(nsg_raw(:,iop));
                   JaveHFES(:,:,iop)    = real(JaveHFES(:,:,iop));
                   JaveHFESPS(:,:,:,iop)    = real(JaveHFESPS(:,:,:,iop));        
               else
                   esgpolylatest_raw(:,:,:,iop) = imag(esgpolylatest_raw(:,:,:,iop));
                   esgpoly_raw(:,:,:,iop) = imag(esgpoly_raw(:,:,:,iop));
                   esgps_raw(:,:,:,iop)   = imag(esgps_raw(:,:,:,iop));
                   esps_raw(:,:,:,iop)    = imag(esps_raw(:,:,:,iop));
                   poly_raw(:,:,iop)      = imag(poly_raw(:,:,iop));
                   ps_raw(:,:,iop)        = imag(ps_raw(:,:,iop));
                   es_raw(:,:,iop)        = imag(es_raw(:,:,iop));
                   esg_raw(:,:,iop)       = imag(esg_raw(:,:,iop));
                   nsg_raw(:,iop)         = imag(nsg_raw(:,iop));
                   JaveHFES(:,:,iop)      = imag(JaveHFES(:,:,iop));
                   JaveHFESPS(:,:,:,iop)  = imag(JaveHFESPS(:,:,:,iop));


              end
          end



  % Avg over rhs's
  for iop=1:5
   % nsg_raw_AVG(iruns,iop) = mean(nsg_raw(:,iop))/hvol;
 nsgAVG(iop) = mean(nsg_raw(:,iop))/hvol;
 nsgERR(iop) = std(nsg_raw(:,iop))/(hvol*sqrt(nrhs));

    for ieig=1:neigH
     % esg_raw_AVG(iruns,ieig,iop) = mean(esg_raw(:,ieig,iop))/hvol;
      esgAVG(ieig,iop) = mean(esg_raw(:,ieig,iop))/hvol;
  esgERR(ieig,iop) = std(esg_raw(:,ieig,iop))/(hvol*sqrt(nrhs));
 
 for ksub=1:nsub+1
      %  esgps_raw_AVG(iruns,ieig,ksub,iop) = mean(esgps_raw(:,ieig,ksub,iop)) / hvol;
 esgpsAVG(ieig,ksub,iop) = mean(esgps_raw(:,ieig,ksub,iop)) / hvol;
 esgpsERR(ieig,ksub,iop) = std(esgps_raw(:,ieig,ksub,iop)) / (hvol*sqrt(nrhs));
JaveHFESPS(ieig,ksub,iop)=std(JaveHFESPS(:,ieig,ksub,iop)) / (hvol*sqrt(nrhs));
% This part is added for Hermitian forced polynomial subtraction BS 5/6/2015
    
esgpolylatestAVG(ieig,ksub,iop) = mean(esgpolylatest_raw(:,ieig,ksub,iop)) / hvol;
esgpolyAVG(ieig,ksub,iop) = mean(esgpoly_raw(:,ieig,ksub,iop)) / hvol;
 esgpolylatestERR(ieig,ksub,iop) = std(esgpolylatest_raw(:,ieig,ksub,iop)) / (hvol*sqrt(nrhs));  
 esgpolyERR(ieig,ksub,iop) = std(esgpoly_raw(:,ieig,ksub,iop)) / (hvol*sqrt(nrhs));  


end
    end

    for ieig=1:neig
 %     es_raw_AVG(iruns,ieig,iop) = mean(es_raw(:,ieig,iop))/hvol;
 esAVG(ieig,iop) = mean(es_raw(:,ieig,iop))/hvol;
 esERR(ieig,iop) = std(es_raw(:,ieig,iop))/(hvol*sqrt(nrhs));
      for ksub=1:nsub+1
       % esps_raw_AVG(iruns,ieig,ksub,iop) = mean(esps_raw(:,ieig,ksub,iop)) / hvol;
 espsAVG(ieig,ksub,iop) = mean(esps_raw(:,ieig,ksub,iop)) / hvol;
 espsERR(ieig,ksub,iop) = std(esps_raw(:,ieig,ksub,iop)) /(hvol*sqrt(nrhs));
     
 end
    end

    for ksub=1:nsub+1
  %    ps_raw_AVG(iruns,ksub,iop) = mean(ps_raw(:,ksub,iop))/hvol;
 psAVG(ksub,iop) = mean(ps_raw(:,ksub,iop))/hvol;
       psERR(ksub,iop) = std(ps_raw(:,ksub,iop))/(hvol*sqrt(nrhs));
    end

% Added part for poly BS 4/3/2015
 for ksub=1:nsub+1
  %    ps_raw_AVG(iruns,ksub,iop) = mean(ps_raw(:,ksub,iop))/hvol;
 polyAVG(ksub,iop) = mean(poly_raw(:,ksub,iop))/hvol;
       polyERR(ksub,iop) = std(poly_raw(:,ksub,iop))/(hvol*sqrt(nrhs));
    end



  end
end

% Average over all trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for iop=1:5
%  if (iop==5) % scalar (use real)
%    nsgAVG(:,iop) = ones(neig,1) * mean(nsg_raw_AVG(:,iop));
%    nsgERR(:,iop) = ones(neig,1) * std(real(nsg_raw_AVG(:,iop)))/sqrt(ntrials);

 %   for ieig=1:neigH
  %    esgAVG(ieig,iop) = mean(esg_raw_AVG(:,ieig,iop));
   %   esgERR(ieig,iop) = std(real(esg_raw_AVG(:,ieig,iop)))/sqrt(ntrials);
    %  for ksub=1:nsub+1
     %   esgpsAVG(ieig,ksub,iop) = mean(esgpnns_raw_AVG(:,ieig,ksub,iop));
      %  esgpsERR(ieig,ksub,iop) = std(real(esgps_raw_AVG(:,ieig,ksub,iop)))/sqrt(ntrials);
      %end
   % end

   % for ieig=1:neig
   %   esAVG(ieig,iop) = mean(es_raw_AVG(:,ieig,iop));
    %  esERR(ieig,iop) = std(real(es_raw_AVG(:,ieig,iop)))/sqrt(ntrials);
     % for ksub=1:nsub+1
      %  espsAVG(ieig,ksub,iop) = mean(esps_raw_AVG(:,ieig,ksub,iop));
       % espsERR(ieig,ksub,iop) = std(real(esps_raw_AVG(:,ieig,ksub,iop)))/sqrt(ntrials);
      %end
   % end

   % for ksub=1:nsub+1
    %  psAVG(ksub,iop) = mean(ps_raw_AVG(:,ksub,iop));
     % psERR(ksub,iop) = std(real(ps_raw_AVG(:,ksub,iop)))/sqrt(ntrials);
   % end
 % else % gamma_mu (use imag)
  %  nsgAVG(:,iop) = ones(neig,1) * mean(nsg_raw_AVG(:,iop));
   % nsgERR(:,iop) = ones(neig,1) * std(imag(nsg_raw_AVG(:,iop)))/sqrt(ntrials);

   % for ieig=1:neigH
    %  esgAVG(ieig,iop) = mean(esg_raw_AVG(:,ieig,iop));
     % esgERR(ieig,iop) = std(imag(esg_raw_AVG(:,ieig,iop)))/sqrt(ntrials);
      %for ksub=1:nsub+1
       % esgpsAVG(ieig,ksub,iop) = mean(esgps_raw_AVG(:,ieig,ksub,iop));
        %esgpsERR(ieig,ksub,iop) = std(imag(esgps_raw_AVG(:,ieig,ksub,iop)))/sqrt(ntrials);
     % end
    %end

    %for ieig=1:neig
     % esAVG(ieig,iop) = mean(es_raw_AVG(:,ieig,iop));
      %esERR(ieig,iop) = std(imag(es_raw_AVG(:,ieig,iop)))/sqrt(ntrials);
      %for ksub=1:nsub+1
       % espsAVG(ieig,ksub,iop) = mean(esps_raw_AVG(:,ieig,ksub,iop));
        %espsERR(ieig,ksub,iop) = std(imag(esps_raw_AVG(:,ieig,ksub,iop)))/sqrt(ntrials);
      %end
    %end

   % for ksub=1:nsub+1
    %  psAVG(ksub,iop) = mean(ps_raw_AVG(:,ksub,iop));
    %  psERR(ksub,iop) = std(imag(ps_raw_AVG(:,ksub,iop)))/sqrt(ntrials);
    %end
 % end
%end
% Store all noise/solution vectors
z=znoise; gx=gxsol;

clear znoise gxsol xsol;
clear iruns rhs ksub xp xpp ieig iop;
clear Minvpert_z cyclemax eig dinv dk ;
clear gdinv gdk gvk invd_vrdag_z invgd_gvrdag_z;

%clear znoise gxsol xsol;
%clear iruns rhs ksub xp xpp ieig iop;
%clear Minvpert_z cyclemax eig dinv dk es_raw es_raw_AVG;
%clear esps_raw_AVG esgps_raw_AVG;
%clear esg_raw esg_raw_AVG esps_raw esgps_raw;
%clear gdinv gdk gvk invd_vrdag_z invgd_gvrdag_z;

clear kB kDR mDR ;
%clear kB kDR mDR ns_raw ns_raw_AVG ps_raw ps_raw_AVG;
