%cd /home/qcd2/matlabtests/hermitian/trace/;
%matlabpool close force local
%matlabpool
%parpool(2)

  latsize = 8888;   % lattice dimensions (string)
  k = 0.1600;       % kappa value
  ntrials = 1;     % trials
  nrhs = 200;      % noises per trial % Chris changed from 200 -> 3
  bc = '1111';      % boundary conditions
  % dimension of quark matrix
  if (latsize == 4444)
    rank = 3072;
  elseif(latsize == 8888)
    rank = 49152;
  end
  nsub = 7;         % maximum subtraction level in PE
  hvol = rank/12;   % hyper volume of lattice

% Chris -- loops over configUSED with cc

  configUSED = 7;   % config used

  %%%%%% DISPLAY RUN INFORMATION %%%%%%%
  disp(['Lattice Size:        ',num2str(latsize)]);
  disp(['Configuration:       ',num2str(configUSED)]);
  disp(['Boundary Conditions: ',bc]);
  disp(['Rank:                ',num2str(rank)]);
  disp(['k:                   ',num2str(k)]);
  disp('');
  disp(['Number of Trials:    ',num2str(ntrials)]);
  disp(['Number of Z4 noise:  ',num2str(nrhs)]);
  disp(['Seed:                ',num2str(1256+(configUSED-1)*rank*nrhs*ntrials+1)]);
  disp(['Working Directory:   ',pwd]);
  disp(['Storage Directory:   ','/data/whytet/Matlab/',num2str(latsize),'/bc',bc,'/',num2str(k*10000)]);

  % file extension for config; {'01','02',...}
  if (configUSED < 10)
    ext = ['0',num2str(configUSED)];
  else
    ext = num2str(configUSED);
  end
% making directory BS 5/29.014
mkdir(['/data/whytet/Matlab/',num2str(latsize),'/bc',bc,'/',num2str(k*10000)])

  % move to proper director for final storage
  % rand('seed',1256+(configUSED-1)*rank*nrhs*ntrials+1); % Chris
  
	% Chris -- creates random noise for each cc
	% Chris -- end noise making
load(['/home/whytet/Matlab/z2--7.mat']);

% Chris -- edit file openers/loaders/ for corect directory and make them change with cc
load(['/home/whytet/Matlab/h07.mat']);
latsize
completetrace
  cd(['/data/whytet/Matlab/',num2str(latsize),'/bc',bc,'/',num2str(k*10000)])
  save(['fig',ext,'-',num2str(latsize),'.mat'],'nsgAVG','nsgERR','esAVG','esERR','esgAVG',...
                                               'esgERR','esgpolyAVG','esgpolyERR','esgpolylatestERR','esgpolylatestAVG','polyERR','polyAVG','psAVG','psERR',...
                                               'espsAVG','espsERR','esgpsAVG','esgpsERR','neig','neigH','JaveHFES','JaveHFESPS','invgxi','invgxipoly')
  save(['vecmodes',ext,'-',num2str(latsize),'.mat'],'z','gx','gd','gvr','M','vr','vl','d');

  save(['heigvals.mat'],'heig');

  save(['evectors.mat'],'veig');


%matlabpool close force local
%parpool close force local
%p = gcp;
%delete(gcp('nocreate'))
