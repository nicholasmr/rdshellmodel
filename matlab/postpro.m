% Nicholas Mossor Rathmann <nicholas.rathmann@gmail.com>, 2014-2016

function postpro(model,nsh,lambda,dtype,ftype,fsh, dumppath)

global TASK;

TASK.CALC_PSEUDOS      = 1;
TASK.KEEP_INVAR_SERIES = 0;
TASK.KEEP_FLUX_SERIES  = 0;
TASK.KEEP_Up_Un        = 0;
TASK.PRINT_INVAR_STATS = 0;
TASK.USE_AGGREGATES    = 1;

%-----------------------------------------
% LOAD SIMULATION DATA
%-----------------------------------------

global Up Un;
global m D S; % Return structures
m = struct();

FILE = sprintf('/spare250gb/dumps/%s/dump__m%i_nsh%i_lam%.1f_dtype%i_ftype%i_fsh%i.nc', dumppath,model,nsh,lambda,dtype,ftype,fsh);
for f = {'model','lambda','parity', 'visc', 'viscinv', 'nt','nti','dt','nshtop','nshbot','k', 'dtype','ftype','fsh', 'p_list', 'q_list', 'G','g','eps','xi'}, m.(char(f)) = ncreadatt(FILE,'/', char(f)); end

m.helmod = m.model; 
m.sigmamodel = m.helmod-10*(m.helmod > 10); 
m.mixed = (m.sigmamodel == 5);

m.rng    = m.nshbot:m.nshtop;
m.nsh    = length(m.rng);
m.frcshell = m.fsh;
m.k      = m.k(m.rng);
m.parity = double(m.parity);
m.q_max  = double(m.nshbot-1);
if sum(m.G(:)>0)>1, m.mixedshapes = 1; else m.mixedshapes = 0; end
m.helicalfrc = (ftype < 0);
m.frclow = (m.frcshell < m.nsh*1/3);
m.frcmid = (m.frcshell >= m.nsh*1/3) && (m.frcshell < m.nsh*2/3);
m.frctop = (m.frcshell >= m.nsh*2/3);

% Restore flattened weight structures
m.num_triads = length(m.p_list);
for fld = {'eps', 'xi', 'g'}, m.(char(fld)) = reshape(m.(char(fld)), m.num_triads,4); end

t_0   = 1;
t_end = 2;
Up = complex(ncread(FILE,'up_re',[m.nshbot t_0],[m.nsh t_end]), ncread(FILE,'up_im',[m.nshbot t_0],[m.nsh t_end]))';
Un = complex(ncread(FILE,'un_re',[m.nshbot t_0],[m.nsh t_end]), ncread(FILE,'un_im',[m.nshbot t_0],[m.nsh t_end]))';
[m.ulen, ~] = size(Up);
m.T   = 1:m.ulen; 
m.NSH = 1:m.nsh;
m.totalsteps = double(m.nt)*double(m.nti);
m.aggregated_steps = double(ncreadatt(FILE,'/', 'AGGREGATION_STEPS')) * double(m.nti);

% Initial spectrum
tref = 1;
m.Up_0 = complex(ncread(FILE,'up_re',[m.nshbot tref],[m.nsh tref]), ncread(FILE,'up_im',[m.nshbot tref],[m.nsh tref]))';
m.Un_0 = complex(ncread(FILE,'un_re',[m.nshbot tref],[m.nsh tref]), ncread(FILE,'un_im',[m.nshbot tref],[m.nsh tref]))';
m.E1p_0 = m.Up_0.*conj(m.Up_0); 
m.E1n_0 = m.Up_0.*conj(m.Up_0);
m.E1_0  = m.E1p_0 + m.E1n_0;

fprintf('----------------------------\n')
fprintf('postpro(): Loaded simulation %s:\n',FILE);
fprintf('\t Model %i \n\t lambda %.2f \n\t nu     %.2e \n\t nu_L   %.2e \n\t FTYPE  %i \n\t DTYPE  %i \n\t sampled tsteps %i \n\t time range load [t_i; t_end] = [t_%i; t_%i]\n', m.sigmamodel, m.lambda, m.visc, m.viscinv, m.ftype, m.dtype, m.ulen, t_0, t_end);
fprintf('----------------------------\n')

%-------------------------------------------
% MODEL GENERATOR SET
%-------------------------------------------

D = struct(); % All data
S = [];
L = m.lambda; % Shortcut
% 
% if TESTNEW
%     D.Eflux_sum = ncread(FILE,'Eflux_sum',[1 1 1 1],[inf 4 inf inf]);
%     aggre_steps = double(ncreadatt(FILE,'/', 'AGGREGATION_STEPS')) * double(m.nti);
%     %aggre_steps = 1000000000
%     m.meaneps         = 2*(1+(m.ftype>0));
%     D.Eflux_sum = D.Eflux_sum/(aggre_steps*m.meaneps); 
%     m.totalsteps = aggre_steps;
%     
%     D.structfuncs = ncread(FILE,'structfuncs_sum',[1 1 1],[inf 4, inf])/aggre_steps;
% end

%--------------------
% E type generators
%--------------------
GEN_E(1,[1 2]) = [1 -L];
GEN_E(2,[1 2]) = [1  L*(L-1)/(L+1)];
GEN_E(3,[1 2]) = [1 -L*(L+1)/(L-1)];
GEN_E(4,[1 2]) = [1  L];
GEN_E(5,[1 2]) = [1  1]; % No E2

% For submodel 2 always pick the pseudo-E generator to be the real non-zero solution.
if m.num_triads
    p = m.p_list(1); q = m.q_list(1);
    c = zeros(q,1); c(1) = 1; c(p+1) = m.eps(1,2); c(q+1) = -(c(p+1)+1);
    c = flipud(c);
    r = roots(c);
    a = log(roots(c))/log(L);
    alphareal_idx = find(~(imag(a(:))) & (abs(real(a(:))) > 1e-5) );
    if length(alphareal_idx) ~= 1
        a(alphareal_idx)
        error('Multiple real non-zero submodel 2 generators?');
    end
    GEN_E(2,2) = r(alphareal_idx);
    fprintf('postpro(): Picked alpha=%.2f instead (real non-zero solution)\n',a(alphareal_idx));
    m.alpha = a(alphareal_idx);
end

ALPHA_E = log(GEN_E)/log(L);

%--------------------
% H type generators
%--------------------
GEN_H(1,[1 2]) = [L -1];
GEN_H(2,[1 2]) = [L -(L-1)/(L+1)];
GEN_H(3,[1 2]) = [L  (L+1)/(L-1)];
GEN_H(4,[1 2]) = [L  1];
GEN_H(5,[1 2]) = [L  L]; % No H2
BETA_H = log(GEN_H)/log(L);

% This model's generator set
MYGENR = struct('E1',   GEN_E(m.sigmamodel,1), 'E2',   GEN_E(m.sigmamodel,2), 'H1',  GEN_H(m.sigmamodel,1), 'H2',  GEN_H(m.sigmamodel,2));
MYEXP  = struct('E1', ALPHA_E(m.sigmamodel,1), 'E2', ALPHA_E(m.sigmamodel,2), 'H1', BETA_H(m.sigmamodel,1), 'H2', BETA_H(m.sigmamodel,2));
fprintf('postpro(): GENERATORS (lambda^alpha_i, lambda^beta_i): E1: %.2f | E2: %.2f | H1: %.2f | H2: %.2f \n', MYGENR.E1, MYGENR.E2, MYGENR.H1, MYGENR.H2);
m.generators = MYGENR;
m.ALPHA_E = ALPHA_E;
m.BETA_H  = BETA_H;

kfrc = m.k(double(m.fsh));
m.meaneps         = 2*(1+(m.ftype>0));
m.meanepspseudo   = kfrc^(MYEXP.E2) * m.meaneps;
m.meandelta       = kfrc^(MYEXP.H1) * m.meaneps;
m.meandeltapseudo = kfrc^(MYEXP.H2) * m.meaneps;

%-------------------------------------------
% SPECTRA AND TIME SERIES OF INVARIANTS
%-------------------------------------------

fprintf('postpro(): Constructing inviscid invariants\n');

D.E1p = ncread(FILE,'up_abs_sum',[m.nshbot],[m.nsh])' ./ m.totalsteps;
D.E1n = ncread(FILE,'un_abs_sum',[m.nshbot],[m.nsh])' ./ m.totalsteps;
D.E1  = D.E1p + D.E1n; 

if TASK.CALC_PSEUDOS
D.E2p = bsxfun(@times, m.k.^(MYEXP.E2), D.E1p);
D.E2n = bsxfun(@times, m.k.^(MYEXP.E2), D.E1n);
D.E2  = D.E2p + D.E2n; D.E2p = []; D.E2n = [];
end

D.H1p = bsxfun(@times, m.k.^(MYEXP.H1), D.E1p);
D.H1n = bsxfun(@times, m.k.^(MYEXP.H1), D.E1n);
D.H1  = D.H1p - D.H1n; D.H1p = []; D.H1n = [];

if TASK.CALC_PSEUDOS
D.H2p = bsxfun(@times, m.k.^(MYEXP.H2), D.E1p);
D.H2n = bsxfun(@times, m.k.^(MYEXP.H2), D.E1n);
D.H2  = D.H2p - D.H2n; D.H2p = []; D.H2n = [];
end

% Sum over all shell contributions (e.g.: E1 = sum_n E_n)
if TASK.KEEP_INVAR_SERIES
    %for INV = {'E1','E2','H1','H2'}
    %for INV = {'E1','H1'}
    for INV = {'E1'}
        INV = char(INV);
        for HEL = {'p','n',''}
            f = strcat(INV,char(HEL));
            D.(strcat(f,'_t')) = sum(D.(f),2); 
        end
    end
end

% For debugging
if TASK.PRINT_INVAR_STATS
    printInviscidStats('E1', 0); 
    printInviscidStats('H1', 1); 
    %printInviscidStats('E2', 0); 
    %printInviscidStats('H2', 1);
end

D.E1_series  = 0; 
D.H1_series  = 0; 
if TASK.KEEP_INVAR_SERIES
    D.E1_series  = D.E1; 
    D.H1_series  = D.H1; 
end

if 1 % ~TASK.USE_AGGREGATES
    D.E1  = mean(D.E1,1)'; 
    D.H1  = mean(D.H1,1)'; 
    if TASK.CALC_PSEUDOS
    D.E2  = mean(D.E2,1)'; 
    D.H2  = mean(D.H2,1)';
    end
end

D.E1p = []; D.E1n = []; 
%D.E1p  = mean(D.E1p,1)'; D.E1n  = mean(D.E1n,1)';

D.Up_end = Up(end,:);
D.Un_end = Un(end,:);
if ~TASK.KEEP_Up_Un, clearvars -global Up Un; end

%-----------------------------------------
% SPECTRAL FLUXES OF INVARIANTS
%-----------------------------------------

D.E1_flux_avg_comp__in_model = ncread(FILE,'Eflux_sum'); 
D.E1_flux_avg_comp__in_model = D.E1_flux_avg_comp__in_model(m.rng,:,:) / (m.meaneps*m.aggregated_steps);

D.E1_flux_avg = zeros(m.nsh);
D.E2_flux_avg = D.E1_flux_avg*0; D.H1_flux_avg = D.E1_flux_avg*0; D.H2_flux_avg = D.E1_flux_avg*0;
PiE1 = zeros(m.nsh,m.num_triads,5); PiE2 = PiE1; PiH1 = PiE1; PiH2 = PiE1;

if TASK.USE_AGGREGATES

    D.corr_p = 2*ncread(FILE,'corr_p_sum'); D.corr_p = D.corr_p(m.rng,:,:);
    D.corr_n = 2*ncread(FILE,'corr_n_sum'); D.corr_n = D.corr_n(m.rng,:,:);

    for ii = 1:m.num_triads
        
        p = m.p_list(ii);
        q = m.q_list(ii);

        CORR_P = squeeze(D.corr_p(:,ii,:) + D.corr_n(:,ii,:));
        CORR_N = squeeze(D.corr_p(:,ii,:) - D.corr_n(:,ii,:));
        
        for jj = [1:4] % Models

            Sp = 1; if (jj == 1) || (jj == 2), Sp = -1; end
            Sq = 1; if (jj == 2) || (jj == 3), Sq = -1; end

            g   = m.g(ii,jj); 
            G   = m.G(ii);
            eps = m.eps(ii,jj);
            if ~(abs(g) > 0), continue; end
            
            %G,g,eps
            %G = 1;
            
            for n = 1:m.nsh

                for mm = [(n+1):(n+q)]

                    if (mm > m.nsh) || ((mm-q) < 1), continue; end % Outside shell range

                    PiE1(n,ii,jj) = PiE1(n,ii,jj) + G*g*m.k(mm-q)^(MYEXP.E1) * m.k(mm-q)*CORR_N(mm-1,jj);
                    PiE2(n,ii,jj) = PiE2(n,ii,jj) + G*g*m.k(mm-q)^(MYEXP.E2) * m.k(mm-q)*CORR_N(mm-1,jj);
                    PiH1(n,ii,jj) = PiH1(n,ii,jj) + G*g*m.k(mm-q)^(MYEXP.H1) * m.k(mm-q)*CORR_P(mm-1,jj);
                    
                    if mm > (n+q-p), continue; end % Second part of sum is over only a sub-range.
                    
                    PiE1(n,ii,jj) = PiE1(n,ii,jj) - Sp*eps* G*g*m.k(mm-q+p)^(MYEXP.E1)* m.k(mm-q)*CORR_N(mm-1,jj);
                    PiE2(n,ii,jj) = PiE2(n,ii,jj) - Sp*eps* G*g*m.k(mm-q+p)^(MYEXP.E2)* m.k(mm-q)*CORR_N(mm-1,jj);
                    PiH1(n,ii,jj) = PiH1(n,ii,jj) - 	 eps* G*g*m.k(mm-q+p)^(MYEXP.H1)* m.k(mm-q)*CORR_P(mm-1,jj);
                end

            end
        end
        
    end
    
    D.E1_flux_avg_comp = PiE1 / (m.meaneps        *m.aggregated_steps);
    D.E2_flux_avg_comp = PiE2 / (m.meanepspseudo  *m.aggregated_steps);
    D.H1_flux_avg_comp = PiH1 / (m.meandelta      *m.aggregated_steps);
    D.H2_flux_avg_comp = PiH2 / (m.meandeltapseudo*m.aggregated_steps);  
    D.E1_flux_avg = sum(sum(D.E1_flux_avg_comp,2),3);
    D.E2_flux_avg = sum(sum(D.E2_flux_avg_comp,2),3);
    D.H1_flux_avg = sum(sum(D.H1_flux_avg_comp,2),3);
end

m.MYEXP=MYEXP;

end

function printInviscidStats(INV, makeabs)

    global D;

    MAXERR = @(Y,Y0) max(abs(Y0-Y));
    RMSE   = @(Y,Y0) sqrt(mean((Y-Y0).^2));
    ABSERR = @(Y,Y0) sum(abs(Y-Y0));
    
    if (makeabs == 1)
        DELIM = '|';
        ABS_ = @(Y) abs(Y);
    else
        DELIM = ' ';
        ABS_ = @(Y) Y;
    end
    
    f = sprintf('%s_t',INV);
    fprintf('------------------------------\n');
    fprintf('%s%s(t)%s MEAN         = %.3e\n', DELIM, INV, DELIM, mean(     ABS_(D.(f))               )  );
    fprintf('%s%s(t)%s MAX ERR      = %.3e\n', DELIM, INV, DELIM, MAXERR(   ABS_(D.(f)),ABS_(D.(f)(1)))  );
    fprintf('%s%s(t)%s RMSE         = %.3e\n', DELIM, INV, DELIM, RMSE(     ABS_(D.(f)),ABS_(D.(f)(1)))  );
    fprintf('%s%s(t)%s SUM ABS ERR  = %.3e\n', DELIM, INV, DELIM, ABSERR(   ABS_(D.(f)),ABS_(D.(f)(1)))  );

   
end