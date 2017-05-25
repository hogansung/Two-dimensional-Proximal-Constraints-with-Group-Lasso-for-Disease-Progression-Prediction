function [] = run_Least_mtl()

addpath(genpath('/Users/hogan/Google Drive/Courses/PreGraduate_1/MALSAR1.1'));

%% set settings
FD = '../dat/TGL_mtl/Longitudinal/';

xList = {'MRI', 'META_MRI'};
yList = {'MMSCORE', 'TOTAL11'};
zList = {'CC', 'RMSE'};

%{
%Least_Lasso_mtl
rLst1(:,:,1) = [0.1, 1; 10, 10];
rLst1(:,:,2) = [0.001, 0.001; 10, 0.001];
rLst2(:,:,1) = [1000, 1000; 1000, 1000];
rLst2(:,:,2) = [1000, 1000; 1000, 1000];
rLst3(:,:,1) = [0.004, 0.004; 0.032, 0.032];
rLst3(:,:,2) = [0.004, 0.004; 0.032, 0.016];
%}

%{
%Least_TGL_mtl_n
rLst1(:,:,1) = [1, 1; 0.001, 0.1];
rLst1(:,:,2) = [0.1, 0.1; 0.001, 0.001];
rLst2(:,:,1) = [1000, 1000; 1000, 1000];
rLst2(:,:,2) = [1000, 1000; 1000, 1000];
rLst3(:,:,1) = [0.001, 0.002; 0.004, 0.008];
rLst3(:,:,2) = [0.001, 0.002; 0.008, 0.008];
%}

%{
%Least_CFGLasso_mtl_n
rLst1(:,:,1) = [0.01, 1; 0.01, 0.1];
rLst1(:,:,2) = [0.01, 1; 0.01, 0.001];
rLst2(:,:,1) = [1000, 1000; 1000, 1000];
rLst2(:,:,2) = [1000, 100; 100, 100];
rLst3(:,:,1) = [0.001, 0.001; 0.004, 0.008];
rLst3(:,:,2) = [0.001, 0.001; 0.004, 0.008];
%}

%{
%Least_TGL_mtl_g
rLst1(:,:,1) = [0.001, 0.001; 0.01, 10];
rLst1(:,:,2) = [1, 0.001; 0.001, 1];
rLst2(:,:,1) = [1000, 1000; 1000, 1000];
rLst2(:,:,2) = [1000, 1000; 1000, 1000];
rLst3(:,:,1) = [0.001, 0.001; 0.001, 0.001];
rLst3(:,:,2) = [0.001, 0.001; 0.001, 0.002];
%}


%Least_CFGLasso_mtl_g_fl
rLst1(:,:,1) = [0.1, 0.001; 0.001, 0.1];
rLst1(:,:,2) = [0.001, 0.1; 0.001, 0.1];
rLst2(:,:,1) = [1000, 100; 1000, 100];
rLst2(:,:,2) = [100, 100; 100, 100];
rLst3(:,:,1) = [0.001, 0.001; 0.001, 0.001];
rLst3(:,:,2) = [0.001, 0.001; 0.001, 0.001];


NL = 4;
nc = 5;
ne = 5;

anames = ['sc '; 'm06'; 'm12'; 'm24'; 'm36'; 'm48'];
t = length(anames) - 1;


for zid = 1:length(zList)
    for yid = 1:length(yList)
        for xid = 1:length(xList)
            
            xname = xList{xid};
            yname = yList{yid};
            zname = zList{zid};
            fprintf('%s %s %s\n', xname, yname, zname);

            %% read data
            cX = cell(t, t);
            cY = cell(t, t);
            for i = 1:t
                for j = i:t%min(i+NL,t)
                    fname = anames(i,:);
                    tname = anames(j+1,:);
                    cY{i,j} = csvread(strcat(FD, xname, '_', yname, '_', fname, '_', tname, '_y.csv'), 1, 0);
                    cX{i,j} = csvread(strcat(FD, xname, '_', yname, '_', fname, '_', tname, '_x.csv'), 1, 0);
                    cX{i,j} = [zscore(cX{i,j}), ones(size(cX{i,j}, 1), 1)];
                end
            end
            d = size(cX{1,1}, 2);
            
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
                idx = cell(t,t);
                for i = 1:t
                    for j = i:min(i+NL,t)
                        n = size(cX{i,j}, 1);
                        idx{i,j} = crossvalind('Kfold', n, nc);
                    end
                end

                tR = zeros(nc,t,t);
                tS = zeros(nc,t,t);
                nw = zeros(nc,t,t);

                for f = 1:nc
                    aX = cell(t,t);
                    aY = cell(t,t);
                    for i = 1:t
                        for j = i:min(i+NL,t)
                            tti = (idx{i,j} == f);
                            tni = ~tti;
                            aX{i,j} = cX{i,j}(tni,:);
                            aY{i,j} = cY{i,j}(tni,:);
                        end
                    end

                    ATy=zeros(d, t, t);
                    for i = 1:t
                        for j = i:min(i+NL,t)
                            ATy(:,i,j)= aX{i,j}'*aY{i,j};
                        end
                    end

                    lambda_max=0;
                    for k=1:d
                        lambda_max=max(lambda_max,...
                            norm( ATy(k,:), 2) );
                    end
                    rho3_n = rho3 * lambda_max;

                    %[w, ~] = Least_Lasso_mtl(aX, aY, rho1, rho2, rho3_n, NL);
                    %[w, ~] = Least_TGL_mtl_n(aX, aY, rho1, rho2, rho3_n, NL);
                    %[w, ~] = Least_CFGLasso_mtl_n(aX, aY, rho1, rho2, rho3_n, NL);
                    %[w, ~] = Least_TGL_mtl_g(aX, aY, rho1, rho2, rho3_n, NL);
                    [w, ~] = Least_CFGLasso_mtl_g_fl(aX, aY, rho1, rho2, rho3_n, NL);

                    for i = 1:t
                        for j = i:min(i+NL,t)
                            tti = (idx{i,j} == f);
                            ttx = cX{i,j}(tti,:);
                            tty = cY{i,j}(tti,:);

                            py = ttx * w(:,i,j);
                            coef = corrcoef(py, tty);
                            tR(f,i,j) = coef(1, 2);
                            tS(f,i,j) = sqrt(mean((tty - py).^2));
                            nw(f,i,j) = size(tti,1);
                        end
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
