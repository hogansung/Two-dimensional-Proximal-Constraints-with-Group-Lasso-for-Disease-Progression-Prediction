function [] = run_Least_sgl()

addpath(genpath('/Users/hogan/Google Drive/Courses/PreGraduate_1/MALSAR1.1'));

%% set settings
FD = '../dat/TGL/Longitudinal/';

xList = {'MRI', 'META_MRI'};
yList = {'MMSCORE', 'TOTAL11'};
zList = {'CC', 'RMSE'};

%{ 
% Least_Lasso_sgl
rLst1(:,:,1) = [1, 10; 1, 1];
rLst1(:,:,2) = [0.01, 0.1; 0.001, 0.001];
rLst2(:,:,1) = [1, 1; 1, 1];
rLst2(:,:,2) = [1, 1; 1, 1];
rLst3(:,:,1) = [0.008, 0.008; 0.032, 0.032];
rLst3(:,:,2) = [0.004, 0.004; 0.016, 0.032];
%}

%{
% Least_TGL_sgl
rLst1(:,:,1) = [1, 0.1; 1, 10];
rLst1(:,:,2) = [1, 0.1; 10, 0.1];
rLst2(:,:,1) = [1000, 1000; 1000, 1000];
rLst2(:,:,2) = [1000, 1000; 1000, 1000];
rLst3(:,:,1) = [0.001, 0.001; 0.001, 0.002];
rLst3(:,:,2) = [0.001, 0.001; 0.001, 0.004];
%}

% Least_CFGLasso_sgl
rLst1(:,:,1) = [0.01, 0.01; 0.001, 1];
rLst1(:,:,2) = [0.1, 0.01; 1, 1];
rLst2(:,:,1) = [1000, 100; 1000, 1000];
rLst2(:,:,2) = [100, 100; 100, 100];
rLst3(:,:,1) = [0.001, 0.001; 0.001, 0.002];
rLst3(:,:,2) = [0.001, 0.002; 0.002, 0.004];

nc = 5;
ne = 5;

tnames = ['m06'; 'm12'; 'm24'; 'm36'; 'm48'];


for zid = 1:length(zList)
    for yid = 1:length(yList)
        for xid = 1:length(xList)
            xname = xList{xid};
            yname = yList{yid};
            zname = zList{zid};
            fprintf('%s %s %s\n', xname, yname, zname);
            
            %% read data
            t = length(tnames);
            cX = cell(t, 1);
            cY = cell(t, 1);
            for i = 1:t
                tname = tnames(i, :);
                cY{i} = csvread(strcat(FD, xname, '_', yname, '_', tname, '_y.csv'), 1, 0);
                cX{i} = csvread(strcat(FD, xname, '_', yname, '_', tname, '_x.csv'), 1, 0);
                cX{i} = [zscore(cX{i}), ones(size(cX{i}, 1), 1)];
            end
            d = size(cX{1}, 2);

            %% set model params
            rho1 = rLst1(yid, xid, zid);
            rho2 = rLst2(yid, xid, zid);
            rho3 = rLst3(yid, xid, zid);
            fprintf('Params: rho1=%f, rho2=%f, rho3=%f\n', rho1, rho2, rho3);
            
            %% set partition seed
            ave = 0;
            var = 0;
            
            for sd = 1:ne
                rng(sd);
            
                %% K-fold partition
                idx = cell(t, 1);
                for i = 1:t
                    n = size(cX{i}, 1);
                    idx{i} = crossvalind('Kfold', n, nc);
                end
          
                tR = zeros(nc,t);
                tS = zeros(nc,t);
                nw = zeros(nc,t);
                for f = 1:nc
                    aX = cell(t,1);
                    aY = cell(t,1);
                    for i = 1:t
                        tti = (idx{i} == f);
                        tni = ~tti;
                        aX{i} = cX{i}(tni,:);
                        aY{i} = cY{i}(tni,:);
                    end

                    ATy=zeros(d, t);
                    for i = 1:t
                        ATy(:,i)= aX{i}'*aY{i};
                    end

                    lambda_max=0;
                    for k=1:d
                        lambda_max=max(lambda_max,...
                            norm( ATy(k,:), 2) );
                    end
                    rho3_n = rho3 * lambda_max;

                    %[w, ~] = Least_Lasso_sgl(aX, aY, rho1, rho2, rho3_n);
                    %[w, ~] = Least_TGL_sgl(aX, aY, rho1, rho2, rho3_n);
                    [w, ~] = Least_CFGLasso_sgl(aX, aY, rho1, rho2, rho3_n);

                    for i = 1:t
                        tti = (idx{i} == f);
                        ttx = cX{i}(tti,:);
                        tty = cY{i}(tti,:);

                        py = ttx * w(:, i);
                        coef = corrcoef(py, tty);
                        tR(f,i) = coef(1, 2);
                        tS(f,i) = sqrt(mean((tty - py).^2));
                        nw(f,i) = size(tti,1);
                    end
                end

                tR = tR .* nw;
                tS = tS .* nw;

                sub_ave = 0;
                if zid == 1
                    for f = 1:nc
                        sub_ave = sub_ave + (sum(tR(f,:)) / sum(nw(f,:)));
                    end
                else
                    for f = 1:nc
                        sub_ave = sub_ave + (sum(tS(f,:)) / sum(nw(f,:)));
                    end
                end
                sub_ave = sub_ave / nc;
                fprintf('%f\n', sub_ave);
                
                ave = ave + sub_ave;
                var = var + sub_ave^2;
            end
            
            ave = ave / ne;
            var = var / ne - ave^2;
            
            fprintf('Total:\n%f\n%f\n\n', ave, var);
        end
    end
end
