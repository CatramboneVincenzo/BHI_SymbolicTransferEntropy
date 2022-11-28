function [STE, Ex, Ey, Len, Surr] = BHI_STE_wShiftStats(microS1_smth, RR, N_labels1, N_labels2, significance, shift)
% BHI_SymbolicTransferEntropy
    [STE, Ex, Ey, Len] = BHI_SymbolicTransferEntropy(microS1_smth, RR, N_labels1, N_labels2);
        
    if significance
        Nperm = 100;
        STEp = zeros(Nperm,2);
        Exp = zeros(Nperm,1);
        Eyp = zeros(Nperm,1);
%         shift = randi([50 200],Nperm,1); % number of samples to perform shift
%         shift_min = 7; % number of samples to perform shift
        for pp = 1:Nperm
            RRperm = [RR(shift(pp):end) RR(1:(shift(pp)-1))];
%             RRperm = [RR((shift_min+pp):end) RR(1:(shift_min+pp-1))];
            microS1_perm = [microS1_smth((end-shift(pp)):end); microS1_smth(1:(end-shift(pp)-1))];
%             microS1_perm = [microS1_smth((end-shift_min-pp):end); microS1_smth(1:(end-shift_min-pp-1))];
%             RRperm = RR(randperm(length(RR)));
%             microS1_perm = microS1_smth(randperm(length(microS1_smth)));
            [STEp(pp,:), Exp(pp), Eyp(pp), ~] = BHI_SymbolicTransferEntropy(microS1_perm, RRperm, N_labels1, N_labels2);
        end
%         pVs = [sum(STE(1)<STEp(:,1)) sum(STE(2)<STEp(:,2)) sum(Ex<Exp) sum(Ey<Eyp)]./Nperm;
%         vec = [(STE-mean(STEp))./std(STEp) (Ex-mean(Exp))/std(Exp) (Ey-mean(Eyp))/std(Eyp)];
%         pVs = 1-cdf('norm',vec,0,1);
        Surr.STE = STEp;
        Surr.Ex = Exp;
        Surr.Ey = Eyp;
    end
end