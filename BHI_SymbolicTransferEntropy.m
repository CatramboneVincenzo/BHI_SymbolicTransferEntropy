function [STE, SEx, SEy, Len] = BHI_SymbolicTransferEntropy(x, y, N_labels1, N_labels2)
% Brain-Heart Iinterplay Symbolic Transfer Entropy with permutation statistics
% output variables:
% STE = Symbolic Transfer Entropy
% SEx = Self Entropy of variable x
% SEy = Self Entropy of variable y
% Len = length of the input series
% Input variables:
% x = vector of N samples (symbolic vector)
% y = vector of N samples (symbolic vector)
% N_labels1 = # of different symbols that variable x might assume
% N_labels2 = # of different symbols that variable y might assume
%
% For methodological details please refer to:
% System-wise Functional Estimation of Directed Brain-Heart Interplay through Microstate Occurrences
% Catrambone Vincenzo, and Valenza Gaetano
% Transaction of Biomedical Engineering, 2023

    if ~isempty(setdiff(1:N_labels1,unique(x)))
        missing = setdiff(1:N_labels1,unique(x));
        x(N_labels1:N_labels1+length(missing)-1) = missing;
    end

    micro = x;
    rr_sym = Nquantization(y,N_labels2);
    Len = length(micro);
    
    if ~isempty(setdiff(1:N_labels2,unique(rr_sym)))
        missing = setdiff(1:N_labels1,unique(x));
        rr_sym(N_labels2:N_labels2+length(missing)-1) = missing;
    end

    micro_k = tuples_conversion(micro, 3);      % x_k
    aa = zeros(N_labels1,N_labels1,N_labels1);
    micro_k_sym = zeros(N_labels1,1);
    for ii = 1:size(micro_k,1)
        aa(micro_k(ii,1),micro_k(ii,2),micro_k(ii,3)) = 1;
        micro_k_sym(ii,1) = find(aa);
        aa(micro_k(ii,1),micro_k(ii,2),micro_k(ii,3)) = 0;
    end

    rr_sym_k = tuples_conversion(rr_sym, 3);
    aa = zeros(N_labels2,N_labels2,N_labels2);
    rr_sym_k_sym = zeros(N_labels2,1);
    for ii = 1:size(rr_sym_k,1)
        aa(rr_sym_k(ii,1),rr_sym_k(ii,2),rr_sym_k(ii,3)) = 1;
        rr_sym_k_sym(ii,1) = find(aa);
        aa(rr_sym_k(ii,1),rr_sym_k(ii,2),rr_sym_k(ii,3)) = 0;
    end

    p_x = accumarray(micro(4:end),1)/size(micro(4:end),1);
    p_x = p_x(:);
    
    p_y = accumarray(rr_sym(4:end),1)/size(rr_sym(4:end),1);
    p_y = p_y(:);
    
    p_xk = accumarray(micro_k_sym,1)/size(micro_k_sym,1);
    p_xk = p_xk(:);
    
    p_yk = accumarray(rr_sym_k,1)/size(rr_sym_k,1);
    p_yk = p_yk(:);

    p_x_xk = accumarray([micro(4:end) micro_k_sym(1:end-1)],1)/(size(micro,1)-3);
    p_y_yk = accumarray([rr_sym(4:end) rr_sym_k_sym(1:end-1)],1)/(size(micro,1)-3);

    p_xk_yk = accumarray([micro_k_sym rr_sym_k_sym],1)/size(micro,1);
    p_yk_xk = accumarray([rr_sym_k_sym micro_k_sym],1)/size(rr_sym_k_sym,1);

    p_x_xk_yk = accumarray([micro(4:end) micro_k_sym(1:end-1) rr_sym_k_sym(1:end-1)],1)/(size(micro,1)-3);
    p_y_yk_xk = accumarray([rr_sym(4:end) rr_sym_k_sym(1:end-1) micro_k_sym(1:end-1)],1)/(size(rr_sym,1)-3);
      

    STE_y_x = 0; 
    STE_x_y = 0;
    SEx = 0;
    SEy = 0;
    for ii = 1:length(micro_k_sym)-1
        p3 = p_x_xk_yk(micro(ii+3),micro_k_sym(ii),rr_sym_k_sym(ii));
        p21 = p_xk_yk(micro_k_sym(ii),rr_sym_k_sym(ii));
        p1 = p_xk(micro_k_sym(ii));
        p22 = p_x_xk(micro(ii+3),micro_k_sym(ii));
        p2 = p_x(micro(ii+3));
        if (p3*log(p3/p21*p1/p22))>=0
            STE_y_x = STE_y_x + p3*log(p3/p21*p1/p22); 
        end
        SEx = SEx + p22*log(p22/(p1*p2));
        
        p3 = p_y_yk_xk(rr_sym(ii+3),rr_sym_k_sym(ii),micro_k_sym(ii));
        p21 = p_yk_xk(rr_sym_k_sym(ii),micro_k_sym(ii));
        p1 = p_yk(rr_sym_k_sym(ii));
        p22 = p_y_yk(rr_sym(ii+3),rr_sym_k_sym(ii));
        p2 = p_y(rr_sym(ii+3));
        if (p3*log(p3/p21*p1/p22))>=0
            STE_x_y = STE_x_y + p3*log(p3/p21*p1/p22);
        end
        SEy = SEy + p22*log(p22/(p1*p2));
    end

    STE = [STE_y_x STE_x_y];
end

