
function R = mymodel_sftf(pars,sfmat,tfmat)

% parse parameters
A = pars(1); % amplitude
mu_sf = pars(2); % center of SF
sigma_sf = pars(3);  % sigma of SF tuning
mu_tf = pars(4);  % center of TF
sigma_tf = pars(5);   % sigma of TF tuning 
Q = pars(6); % -1 to 1, the tan of x/y

% computer terms
tfp = 2.^((Q+1)*( log2(sfmat) - log2(mu_sf) ) + log2(mu_tf));
t1 = exp( -(log2(sfmat)-log2(mu_sf)).^2 / log2(sigma_sf).^2 );
t2 = exp( -(log2(tfmat)-log2(tfp)).^2 / log2(sigma_tf).^2 );

R = A*t1.*t2;

return;

