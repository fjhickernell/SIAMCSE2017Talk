function [K,kvec,k0] = kernelFun(x,whKer,shape,domain,BernPolynX,BernPolynOrder)
[nx,d] = size(x);
if nargin < 4
   domain = [zeros(1,d); ones(1,d)];
   if nargin < 3;
      shape = 1;
      if nargin < 2
         whKer = 'sqExp';
      end
   end
end

if ~strcmp(whKer,'Fourier')
    K = ones(nx);
end
if strcmp(whKer,'sqExp') %all wrong
   kvec = ones(nx,1)*(sqrt(pi)/(2*shape))^d;
   for k = 1:d;
      K = K.*exp(-(shape*bsxfun(@minus,x(:,k),x(:,k)')).^2);
      kvec = kvec.*(erf(shape*x(:,k)) + erf(shape*(1 - x(:,k))));
   end
elseif strcmp(whKer,'Mat1')
   diffdom = diff(domain,1,1);
   shdiffdom = shape*diffdom;
   k0 = prod((- 6 + 4*shdiffdom +exp(-shdiffdom).*(6 + 2*shdiffdom))./shdiffdom.^2);
   kvec = ones(nx,1)*(2^d/prod(shdiffdom));
   for k = 1:d;
      tempa = shape*abs(bsxfun(@minus,x(:,k),x(:,k)'));
      K = K.*exp(-tempa).*(1 + tempa);
      tempb = shape*(x(:,k)-domain(1,k));
      tempc = shape*(domain(2,k) - x(:,k));
      kvec = kvec.*(2 - exp(-tempc).*(1+tempc/2) ...
          - exp(-tempb).*(1+tempb/2));
   end
elseif strcmp(whKer,'Fourier')
    r = BernPolynOrder;
	theta = shape;
	kvec = ones(nx,1);
	k0 = 1.0;
	constMult = -theta^r*(-1)^(r/2)*(2*pi)^r/factorial(r);
    bernX = BernPolynX*constMult;

    cvec = prod((1.0+bernX),2);
    % no need to explicitly create the matrix
    % circulant matrix can be recovered with just one vector
    K = cvec';
end
    
end





