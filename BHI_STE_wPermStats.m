function [STE, SEx, SEy, pVs] = BHI_STE_wPermStats(x, y, N_labels1, N_labels2, significance, Nperm)
% Brain-Heart Iinterplay Symbolic Transfer Entropy with permutation statistics
% output variables:
% STE = Symbolic Transfer Entropy
% SEx = Self Entropy of variable x
% SEy = Self Entropy of variable y
% pVs = p-value associate to STE
% Input variables:
% x = vector of N samples (symbolic vector)
% y = vector of N samples (symbolic vector)
% N_labels1 = # of different symbols that variable x might assume
% N_labels2 = # of different symbols that variable y might assume
% significance = boolean, if true surrogate data analysis is performed and pVs calculated
% Nperm = # of permutation to calculate p-values
%
% For details refer to:
% System-wise Functional Estimation of Directed Brain-Heart Interplay through Microstate Occurrences
% Catrambone Vincenzo, and Valenza Gaetano
% Transaction of Biomedical Engineering, 2023

    [STE, SEx, SEy, ~] = BHI_SymbolicTransferEntropy(x, y, N_labels1, N_labels2);
        
    if significance
        STEp = zeros(Nperm,2);
        for pp = 1:Nperm
            try
            RRperm = y(randperm(length(y)));
            MICROperm = x(randperm(length(x)));

            [Surr, ~, ~, ~] = BHI_SymbolicTransferEntropy(x, RRperm, N_labels1, N_labels2);
            STEp(pp,1) = Surr(1);
            [Surr, ~, ~, ~] = BHI_SymbolicTransferEntropy(MICROperm, y, N_labels1, N_labels2);
            STEp(pp,2) = Surr(2);
            catch
                pp = pp - 1;
            end
        end

        pVs = min(1-(sum(STE>STEp)/Nperm),(sum(STE>STEp)/Nperm))*2;
    end
end