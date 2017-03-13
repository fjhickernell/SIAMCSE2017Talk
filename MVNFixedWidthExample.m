%% Generate Examples of Multivariate Normal Probabilities

gail.InitializeWorkspaceDisplay %clean up 
format long

C = [4 1 1; 0 1 0.5; 0 0 0.25];
Cov = C'*C
mu = 0;
a = [-6 -2 -2];
b = [5 2 1];
nRep = 100;
alpha = 0.1; %top of quantile
absTol = 0;
relTol = 0.0001;
dataFileName = ['MVNProbFixedWidthExampleAllDataA' ...
   int2str(log10(absTol)) 'R' int2str(log10(relTol)) '.mat'];


if exist(dataFileName,'file')
   load(dataFileName)
   MVNProbBestArch = MVNProbBest;
   nRepGoldArch = nRepGold;
   %MVNProbIIDGgArch = MVNProbIIDGg;
   MVNProbSobolGgArch = MVNProbSobolGg;
   MVNProbMLELatticeGgArch = MVNProbMLELatticeGg;
   if exist('MVNProbMLELatticeGn', 'var')
        MVNProbMLELatticeGnArch = MVNProbMLELatticeGg;
   end
end

%% First compute a high accuracy answer
nGold = 2^27;
nRepGold = nRep;
MVNProbBest = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nGold, ...
   'errMeth','n','cubMeth','Sobol','intMeth','Genz');
compGold = true;
if exist(dataFileName,'file')
   if sameProblem(MVNProbBest,MVNProbBestArch) && ...
      nRepGoldArch == nRepGold
      disp('Already have gold standard answer')
      compGold = false;
   end
end
if false
   disp('(Re-)computing gold standard answer')
   muBestvec = zeros(1,nRepGold);
   tic 
   for i = 1:nRepGold
      i
      tic
      muBestvec(1,i) = compProb(MVNProbBest); 
      toc
   end
   toc
   muBest = mean(muBestvec);
end
disp(['mu  = ' num2str(muBest,15) ' +/- ' num2str(2*std(muBestvec),10)])

%% IID sampling
if absTol + relTol >= 0.001
MVNProbIIDGg = multivarGauss('a',a,'b',b,'Cov',Cov, ...
   'errMeth','g','absTol',absTol,'relTol',relTol,...
   'cubMeth','IID','intMeth','Genz');
compIID = true;
if exist(dataFileName,'file')
   if sameProblem(MVNProbIIDGg,MVNProbIIDGgArch) && ...
      all(nRep == nRepArch)
      disp('Already have IID answer')
      compIID = false;
   end
end
if compIID
   tic
   muMVNProbIIDGg = zeros(nRep,1);
   IIDGgNrequired = zeros(nRep,1);
   IIDGgTime = zeros(nRep,1);
   for i = 1:nRep
      gail.TakeNote(i,10)
      [muMVNProbIIDGg(i),out] = compProb(MVNProbIIDGg);
      IIDGgNrequired(i) = out.ntot;
      IIDGgTime(i) = out.time;
   end
   errvecMVNProbIIDGg = abs(muBest - muMVNProbIIDGg);
   errmedMVNProbIIDGg = median(errvecMVNProbIIDGg);
   errtopMVNProbIIDGg = quantile(errvecMVNProbIIDGg,1-alpha);
   errSucceedIIDGg = mean(errvecMVNProbIIDGg <= max(absTol,relTol*abs(muBest)));
   topNIIDGg = quantile(IIDGgNrequired,1-alpha);
   topTimeIIDGg = quantile(IIDGgTime,1-alpha);
   toc
end
end

%% Scrambled Sobol sampling
MVNProbSobolGg = multivarGauss('a',a,'b',b,'Cov',Cov, ...
   'errMeth','g','absTol',absTol,'relTol',relTol,...
   'cubMeth','Sobol','intMeth','Genz');
compSobol = true;
if exist(dataFileName,'file')
   if sameProblem(MVNProbSobolGg,MVNProbSobolGgArch)
      disp('Already have Scrambled Sobol answer')
      compSobol = false;
   end
end
if true
   tic 
   muMVNProbSobolGg = zeros(nRep,1);
   SobolGgNrequired = zeros(nRep,1);
   SobolGgTime = zeros(nRep,1);
   for i = 1:nRep
      [muMVNProbSobolGg(i), out] = compProb(MVNProbSobolGg); 
      SobolGgNrequired(i) = out.n;
      SobolGgTime(i) = out.time;
   end
   errvecMVNProbSobolGg = abs(muBest - muMVNProbSobolGg);
   errmedMVNProbSobolGg = median(errvecMVNProbSobolGg);
   errtopMVNProbSobolGg = quantile(errvecMVNProbSobolGg,1-alpha);
   errSucceedSobolGg = mean(errvecMVNProbSobolGg <= max(absTol,relTol*abs(muBest)));
   topNSobolGg = quantile(SobolGgNrequired,1-alpha);
   topTimeSobolGg = quantile(SobolGgTime,1-alpha);
   toc
end

%% Try MLE Bayseian cubature with Fourier kernel and Rank1 Lattice points
MVNProbMLELatticeGg = multivarGauss('a',a,'b',b,'Cov',Cov, ...
   'errMeth','g','absTol',absTol,'relTol',relTol,...
   'cubMeth','MLELattice','intMeth','Genz');
compMLELattice = true;
if exist(dataFileName,'file') % force to compute all the time
   if sameProblem(MVNProbMLELatticeGg,MVNProbMLELatticeGgArch)
      disp('Already have MLE Fourier Lattice answer')
      compMLELattice = false;
   end
end
if true
   %datetime
   tic
   muMVNProbMLELatticeGg = zeros(nRep,1);
   MLELatticeGgNrequired = zeros(nRep,1);
   MLELatticeGgTime = zeros(nRep,1);
   errbdvecMBVProbMLELatticeGg(nnMLE,nRep) = 0;
   for i = 1:nRep
      gail.TakeNote(i,10)
      [muMVNProbMLELatticeGg(:,i), out] = compProb(MVNProbMLELatticeGg); 
      errbdvecMBVProbMLELatticeGg(:,i) = out.ErrBd;
      MLELatticeGgNrequired(i) = out.n;
      MLELatticeGgTime(i) = out.time;
   end
   errvecMVNProbMLELatticeGg = abs(muBest - muMVNProbMLELatticeGg);
   errmedMVNProbMLELatticeGg = median(errvecMVNProbMLELatticeGg);
   errtopMVNProbMLELatticeGg = quantile(errvecMVNProbMLELatticeGg,1-alpha);
   errSucceedMLELatticeGg = mean(errvecMVNProbSobolGg <= max(absTol,relTol*abs(muBest)));
   topNMLELatticeGg = quantile(MLELatticeGgNrequired,1-alpha);
   topTimeMLELatticeGg = quantile(MLELatticeGgTime,1-alpha);
   toc
   datetime
end

%% Save output
save(dataFileName)
return

