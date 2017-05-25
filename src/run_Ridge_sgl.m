function [] = run_Ridge_sgl()

addpath(genpath('/Users/hogan/Google Drive/Courses/PreGraduate_1/MALSAR1.1'));
warning('off');

%% set settings
FD = '../dat/TGL/Longitudinal/';

xList = {'MRI', 'META_MRI'};
yList = {'MMSCORE', 'TOTAL11'};
zList = {'CC', 'RMSE'};

rLst1(:,:,1) = [1000, 1000; 0.001, 0.001];
rLst1(:,:,2) = [1000, 1000; 0.001, 0.001];

rLst2(:,:,1) = [1, 1; 1, 1];
rLst2(:,:,2) = [1, 1; 1, 1];

rLst3(:,:,1) = [1, 1; 1, 1];
rLst3(:,:,2) = [1, 1; 1, 1];

nc = 5;
ne = 5;

for zid = 1:length(zList)
    for yid = 1:length(yList)
        for xid = 1:length(xList)
            xname = xList{xid};
            yname = yList{yid};
            zname = zList{zid};
            fprintf('%s %s %s\n', xname, yname, zname);
            
            %% read data
            cX = cell(1, 1);
            cY = cell(1, 1);
            cY{1} = csvread(strcat(FD, xname, '_', yname, '_y.csv'), 1, 0);
            cX{1} = csvread(strcat(FD, xname, '_', yname, '_x.csv'), 1, 0);
            cX{1} = [zscore(cX{1}), ones(size(cX{1}, 1), 1)];
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
                idx = cell(1, 1);
                n = size(cX{1}, 1);
                idx{1} = crossvalind('Kfold', n, nc);

                %% main
                tR = zeros(nc,1);
                tS = zeros(nc,1);
                nw = zeros(nc,1);
                for f = 1:nc
                    aX = cell(1,1);
                    aY = cell(1,1);
                    tti = (idx{1} == f);
                    tni = ~tti;
                    aX{1} = cX{1}(tni,:);
                    aY{1} = cY{1}(tni,:);

                    w = (aX{1}'*aX{1} + 2 * rho1) \ (aX{1}'*aY{1});

                    tti = (idx{1} == f);
                    ttx = cX{1}(tti,:);
                    tty = cY{1}(tti,:);

                    py = ttx * w(:, 1);
                    coef = corrcoef(py, tty);
                    tR(f,1) = coef(1, 2);
                    tS(f,1) = sqrt(mean((tty - py).^2));
                    nw(f,1) = size(tti,1);
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
