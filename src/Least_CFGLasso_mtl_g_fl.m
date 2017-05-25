%% FUNCTION Least_CFGLasso
%   Fused Group Lasso Progression Model with Least Squares Loss.
%
%% OBJECTIVE
%   argmin_W { sum_i^t (0.5 * norm (Y{i} - X{i}' * W(:, i))^2)
%            + rho1 * \|W\|_1 + rho2 * \|W*R\|_1  +
%            + rho3 * \|W\|_{2,1} }
%   where
%   rho1: sparse.
%   rho2: fused Lasso.
%   rho3: L2,1-norm.
%   R encodes fused structure relationship [1 -1 0 ...; 0 1 -1 ...; ...]
%      R=zeros(t,t-1);R(1:(t+1):end)=1;R(2:(t+1):end)=-1;
%
%% INPUT
%   X: {n * d} * t - input matrix
%   Y: {n * 1} * t - output matrix
%   rho1: sparse.
%   rho2: fused Lasso.
%   rho3: L2,1-norm.
%
%% OUTPUT
%   W: model: d * t
%   funcVal: function value vector.
%
%% LICENSE
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Copyright (C) 2011 - 2012 Jiayu Zhou and Jieping Ye
%
%   You are suggested to first read the Manual.
%   For any problem, please contact with Jiayu Zhou via jiayu.zhou@asu.edu
%
%   Last modified on June 3, 2012.
%
%% RELATED PAPERS
%
%   [1] Zhou, J., Jun, L., Narayan, A. V. and Ye, J.  Modeling Disease
%   Progression via Fused Sparse Group Lasso. KDD 2012
%
%% RELATED FUNCTIONS
%   Logistic_CFGLasso, init_opts

function [W, funcVal] = Least_CFGLasso_mtl_g_fl(X, Y, rho1, rho2, rho3, NL, opts)

if nargin < 6
    error('\n Inputs: X, Y, rho1, rho2, rho3 and NL should be specified!\n');
end

task_num  = size (X, 1);
for t_idx = 1:task_num
    for t_jdx = t_idx:min(t_idx+NL,task_num)
        X{t_idx,t_jdx} = X{t_idx,t_jdx}';
    end
end

if nargin < 7
    opts = [];
end

% initialize options.
opts=init_opts(opts);

dimension = size(X{1,1}, 1);
funcVal = [];

% Relation
% R=zeros(task_num,task_num-1);
% R(1:(task_num+1):end)=1;
% R(2:(task_num+1):end)=-1;
% R = R';

R=cell(task_num,1);
for t_idx = 1:task_num
    R{t_idx} = zeros(task_num,task_num-1);
    for t_kdx = t_idx:min(t_idx+NL,task_num)-1
        R{t_idx}(t_kdx,t_kdx) = 1;
        R{t_idx}(t_kdx+1,t_kdx) = -1;
    end
end

S=cell(task_num,1);
for t_jdx = 1:task_num
    S{t_jdx} = zeros(task_num,task_num-1);
    for t_kdx = max(t_jdx-NL, 1):t_jdx-1
        S{t_jdx}(t_kdx,t_kdx) = 1;
        S{t_jdx}(t_kdx+1,t_kdx) = -1;
    end
end

% initial W
%W0 = zeros(dimension, task_num);
if opts.init==2
    W0 = zeros(dimension, task_num, task_num);
elseif opts.init== 0
    W0 = randn(dimension, task_num, task_num);
else
    if isfield(opts,'W0')
        W0=opts.W0;
        if (nnz(size(W0)-[dimension, task_num, task_num]))
            error('\n Check the input .W0');
        end
    else
        W0 = zeros(dimension, task_num, task_num);
    end
end

bFlag=0; % this flag tests whether the gradient step only changes a little


Wz= W0;
Wz_old = W0;

t = 1;
t_old = 0;

iter = 0;
gamma = 1;
gamma_inc = 2;

while iter < opts.maxIter
    alpha = (t_old - 1) /t;
    
    Ws = (1 + alpha) * Wz - alpha * Wz_old;
    
    % compute function value and gradients of the search point
    gWs  = gradVal_eval(Ws);
    Fs   = funVal_eval (Ws);
    
    while true
        Wzp = FGLasso_projection(Ws - gWs/gamma, rho1 / gamma, rho2 / gamma, rho3 / gamma);
        Fzp = funVal_eval  (Wzp);
        
        delta_Wzp = Wzp - Ws;
        r_sum = 0;
        for t_idx = 1:task_num
            for t_jdx = t_idx:task_num
                r_sum = r_sum + norm(delta_Wzp(:,t_idx,t_jdx))^2;
            end
        end
        
        Fzp_gamma = Fs + sum(sum(sum(delta_Wzp .* gWs)))...
            + gamma/2 * r_sum;
        
        if (r_sum <= 1e-20)
            bFlag=1; % this shows that, the gradient step makes little improvement
            break;
        end
        
        if (Fzp <= Fzp_gamma)
            break;
        else
            gamma = gamma * gamma_inc;
        end
    end
    
    Wz_old = Wz;
    Wz = Wzp;
    
    funcVal = cat(1, funcVal, Fzp + nonsmooth_eval(Wz, rho1, rho2, rho3));
    
    if (bFlag)
        % fprintf('\n The program terminates as the gradient step changes the solution very small.');
        break;
    end
    
    % test stop condition.
    switch(opts.tFlag)
        case 0
            if iter>=2
                if (abs( funcVal(end) - funcVal(end-1) ) <= opts.tol)
                    break;
                end
            end
        case 1
            if iter>=2
                if (abs( funcVal(end) - funcVal(end-1) ) <=...
                        opts.tol* funcVal(end-1))
                    break;
                end
            end
        case 2
            if ( funcVal(end)<= opts.tol)
                break;
            end
        case 3
            if iter>=opts.maxIter
                break;
            end
    end
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
end

W = Wzp;

% private functions

    function [Wp] = FGLasso_projection (W, lambda_1, lambda_2, lambda_3 )
        % solve it in row wise (L_{2,1} is row coupled).
        % for each row we need to solve the proximal opterator
        % argmin_w { 0.5 \|w - v\|_2^2
        %            + lambda_1 * \|w\|_1 + lambda_2 * \|R * w\|_1  +
        %            + lambda_3 * \|w\|_2 }
        % NOTE: Here the R is t-1 * t, the outside R is t * t-1, and
        
        Wp = zeros(size(W));
        
        for k = 1:size(W, 1)
            v = squeeze(W(k,:,:));
            w = FGLasso_projection_full(v, lambda_1, lambda_2, lambda_3);
            Wp(k,:,:) = w;
        end
    end

    function w = FGLasso_projection_full(v, lambda_1, lambda_2, lambda_3)
        w_1 = TDFusedLasso(v, size(v, 1), lambda_1, lambda_2);
        nm = norm(w_1, 'fro');
        if nm == 0
            w_2 = zeros(size(w_1));
        else
            w_2 = max(nm - lambda_3, 0)/nm * w_1;
        end
        w = w_2;
    end

% smooth part gradient.
    function [grad_W] = gradVal_eval(W)
        grad_W = zeros(size(W));
        for i = 1:task_num
            for j = i:min(i+NL,task_num)
                grad_W(:,i,j) = X{i,j}*(X{i,j}' * W(:,i,j)-Y{i,j});
            end
        end
    end

% smooth part gradient.
    function [funcVal] = funVal_eval (W)
        funcVal = 0;
        for i = 1:task_num
            for j = i:min(i+NL,task_num)
                funcVal = funcVal + 0.5 * norm (Y{i,j} - X{i,j}' * W(:,i,j))^2;
            end
        end
    end

    function [non_smooth_value] = nonsmooth_eval(W, rho_1, rho_2, rho_3)
        non_smooth_value = 0;
        for k = 1 : size(W, 1)
            non_smooth_value = non_smooth_value ...
                + rho_1 * norm(squeeze(W(k,:,:)), 1) ...
                + rho_3 * norm(squeeze(W(k,:,:)), 'fro');
            for i = 1:task_num
                non_smooth_value = non_smooth_value ...
                    + rho_2 * norm(squeeze(W(k,i,:))' * R{i}, 1);
            end
            for j = 1:task_num
                non_smooth_value = non_smooth_value ...
                    + rho_2 * norm(squeeze(W(k,:,j)) * S{j}, 1);
            end
        end
    end

end