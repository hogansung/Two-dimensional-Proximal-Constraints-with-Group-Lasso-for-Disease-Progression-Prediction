FD = '../dat_sgl/TGL/Longitudinal/';

xLst = ['MRI', 'META_MRI']
yLst = ['MMSCORE', 'TOTAL11']
tLst = ['m06', 'm12', 'm24', 'm36', 'm48']
sLst = ['x', 'y']

for xn in xLst:
    for yn in yLst:
        for sn in sLst:
            fnn = FD + '_'.join([xn, yn, sn]) + '.csv'

            with open(fnn, 'w') as otf:
                fst_flag = True

                for tn in tLst:
                    fon = FD + '_'.join([xn, yn, tn, sn]) + '.csv'

                    with open(fon) as inf:
                        lines = inf.readlines()
                        if fst_flag is True:
                            otf.write(''.join(lines))
                            fst_flag = False
                        else:
                            otf.write(''.join(lines[1:]))
