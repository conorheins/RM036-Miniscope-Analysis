function [ A,C,C_raw,S,kernel_pars,smin,sn ] = updateTemporal_endoscope_noObj(Y,A,C,numIters,deconv_flag,deconv_options)
%% run HALS by fixating all spatial components
% input:
%   Y:  d*T, fluorescence data
%   allow_deletion: boolean, allow deletion (default: true) 
% output:
%   C_raw: K*T, temporal components without being deconvolved

% options
maxIter = numIters;
deconv_options_0 = deconv_options;


%% initialization
K = size(A, 2);     % number of components
C_raw = zeros(size(C));
S = zeros(size(C));
A = full(A);
U = A'*Y;
V = A'*A;
aa = diag(V);   % squares of l2 norm all all components
sn =  zeros(1, K);
smin = zeros(1,K);
kernel_pars = cell(K,1);

%% updating
for miter=1:maxIter
    for k=1:K
        temp = C(k, :) + (U(k, :)-V(k, :)*C)/aa(k);
        %remove baseline and estimate noise level
        if range(temp)/std(temp)>6
            [b, tmp_sn] = estimate_baseline_noise(temp);
        else
            b = mean(temp(temp<median(temp)));
            tmp_sn = GetSn(temp);
        end
        % we use two methods for estimating the noise level
        %         psd_sn = GetSn(temp);
        %         if tmp_sn>psd_sn
        %             tmp_sn =psd_sn;
        %             [temp, ~] = remove_baseline(temp, tmp_sn);
        %         else
        %             temp = temp - b;
        %         end
        temp = temp -b;
        sn(k) = tmp_sn;
        
        % deconvolution
        if deconv_flag
            [ck, sk, deconv_options]= deconvolveCa(temp, deconv_options_0, 'maxIter', 2);
            smin(k) = deconv_options.smin;
            kernel_pars{k} = reshape(deconv_options.pars, 1, []);
        else
            ck = max(0, temp);
        end
        % save convolution kernels and deconvolution results
        C(k, :) = ck;
        
        % save the spike count in the last iteration
        if miter==maxIter
            if deconv_flag
                S(k, :) = sk;
            end
            C_raw(k, :) = temp;
        end
    end
end

A = bsxfun(@times, A, sn);
C = bsxfun(@times, C, 1./sn');
C_raw = bsxfun(@times, C_raw, 1./sn');
S = bsxfun(@times, S, 1./sn');
kernel_pars =cell2mat( kernel_pars);
smin = smin/sn;

end

