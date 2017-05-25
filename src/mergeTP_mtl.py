FD = '../dat/TGL_mtl/Longitudinal/';

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

                for i in range(len(tLst)):
                    for j in range(i+1, len(tLst)):
                        tn = tLst[i]
                        tn2 = tLst[j]

                        fon = FD + '_'.join([xn, yn, tn, tn2, sn]) + '.csv'

                        with open(fon) as inf:
                            lines = inf.readlines()
                            if fst_flag is True:
                                otf.write(''.join(lines))
                                fst_flag = False
                            else:
                                otf.write(''.join(lines[1:]))
