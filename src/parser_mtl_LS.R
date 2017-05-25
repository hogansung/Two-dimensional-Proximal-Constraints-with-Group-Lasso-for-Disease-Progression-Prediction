library(plyr)


FD = '../dat/Origin/'
TD = '../dat/TGL_mtl/Longitudinal/'


# read files
mri_all = read.csv(paste(FD, 'UCSFFSL_02_01_16.csv', sep=''))

# mri_all
mri_all = mri_all[, colSums(is.na(mri_all)) != nrow(mri_all)] # remove all NA columns
mri_all = mri_all[mri_all$STATUS == 'complete', ] # only complete mri_all
mri_all = mri_all[mri_all$OVERALLQC == 'Pass', ]
drops = c('VISCODE2','EXAMDATE','VERSION','FLDSTRENG','LONISID','LONIUID','IMAGEUID','RUNDATE','STATUS','BASETP1','BASETP2','BASETP3','BASETP4','BASETP5','BASETP6','BASETP7','BASETP8','OVERALLQC','TEMPQC','FRONTQC','PARQC','INSULAQC','OCCQC','BGQC','CWMQC','VENTQC', 'update_stamp','ST8SV') # remain 'VISCODE'
mri_all = mri_all[, !colnames(mri_all) %in% drops]
mri_all = mri_all[complete.cases(mri_all), ]

# mri 
mri = mri_all[mri_all$VISCODE == 'sc', 'RID', drop=FALSE]

# merge MMSE baseline score
mmse = read.csv(paste(FD, 'MMSE.csv', sep=''))
mmse = mmse[(mmse$Phase=='ADNI1') & (mmse$VISCODE=='sc'), ]
mri_mmse = merge(mri, mmse)
mri_mmse$sc_MMSCORE = mri_mmse$MMSCORE
nnames = c('sc_MMSCORE')
for (nname in nnames) mri_mmse = mri_mmse[!(mri_mmse[[nname]] %in% c(-1,-4)),]
mri = mri_mmse[, c(colnames(mri), nnames)]

# Data: MMSE Features
mri = mri[complete.cases(mri), ]
mri_x = mri
mri_fnames = colnames(mri)[-1] # remove RID

# merge ADAS baseline score
adas = read.csv(paste(FD, 'ADASSCORES.csv', sep=''))
adas = adas[adas$VISCODE=='bl', ]
mri_adas = merge(mri, adas)
onames = c('TOTAL11','TOTALMOD','Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8','Q9','Q10','Q11','Q12','Q14')
nnames = sapply(onames, function(x) paste('sc_', x, sep=''))
for (i in 1:length(onames)) mri_adas[[nnames[i]]] = mri_adas[[onames[i]]]
for (nname in nnames) mri_adas = mri_adas[!(mri_adas[[nname]] %in% c(-1,-4)),]
mri = mri_adas[, c(colnames(mri), nnames)]

# merge Demonstration score
dem = read.csv(paste(FD, 'PTDEMOG.csv', sep=''))
dem = dem[(dem$Phase=='ADNI1') & (dem$VISCODE=='sc'), ]
dem$PTAGE = 2005 - dem$PTDOBYY
mri_dem = merge(mri, dem)
nnames = c('PTAGE', 'PTEDUCAT', 'PTGENDER')
for (nname in nnames) mri_dem = mri_dem[!(mri_dem[[nname]] %in% c(-1,-4)),]
mri = mri_dem[, c(colnames(mri), nnames)]

# merge Clinical Dementia Rating Scale score
cdr = read.csv(paste(FD, 'CDR.csv', sep=''))
cdr = cdr[(cdr$Phase=='ADNI1') & (cdr$VISCODE=='sc'), ]
mri_cdr = merge(mri, cdr)
nnames = c('CDMEMORY','CDORIENT','CDJUDGE','CDCOMMUN','CDHOME','CDCARE','CDGLOBAL')
for (nname in nnames) mri_cdr = mri_cdr[!(mri_cdr[[nname]] %in% c(-1,-4)),]
mri = mri_cdr[, c(colnames(mri), nnames)]

# merge Functional Activities Questionnaire score
faq = read.csv(paste(FD, 'FAQ.csv', sep=''))
faq = faq[(faq$Phase=='ADNI1') & (faq$VISCODE=='bl'), ]
mri_faq = merge(mri, faq)
nnames = c('FAQTOTAL')
for (nname in nnames) mri_faq = mri_faq[!(mri_faq[[nname]] %in% c(-1,-4)),]
mri = mri_faq[, c(colnames(mri), nnames)]

# merge Geriatric Depression Scale score
gds = read.csv(paste(FD, 'GDSCALE.csv', sep=''))
gds = gds[(gds$Phase=='ADNI1') & (gds$VISCODE=='sc'), ]
mri_gds = merge(mri, gds)
nnames = c('GDTOTAL')
for (nname in nnames) mri_gds = mri_gds[!(mri_gds[[nname]] %in% c(-1,-4)),]
mri = mri_gds[, c(colnames(mri), nnames)]

# merge Modified Hachinski Ischemia Scale score
hach = read.csv(paste(FD, 'MODHACH.csv', sep=''))
hach = hach[(hach$Phase=='ADNI1') & (hach$VISCODE=='sc'), ]
mri_hach = merge(mri, hach)
nnames = c('HMSCORE')
for (nname in nnames) mri_hach = mri_hach[!(mri_hach[[nname]] %in% c(-1,-4)),]
mri = mri_hach[, c(colnames(mri), nnames)]

# merge Neuropsychological Battery score
nb = read.csv(paste(FD, 'NEUROBAT.csv', sep=''))
nb = nb[(nb$Phase=='ADNI1') & (nb$VISCODE=='sc'), ]
mri_nb = merge(mri, nb)
nnames = c('LIMMTOTAL')
for (nname in nnames) mri_nb = mri_nb[!(mri_nb[[nname]] %in% c(-1,-4)),]
mri = mri_nb[, c(colnames(mri), nnames)]

# merge Labortary score
lab = read.csv(paste(FD, 'LABDATA.csv', sep=''))
lab = lab[(lab$Phase=='ADNI1') & (lab$VISCODE=='sc'), ]
lab = data.frame(lapply(lab, function(x) suppressWarnings(as.numeric(as.character(x)))))
mri_lab = merge(mri, lab)
nnames = c('RCT1','RCT11','RCT12','RCT13','RCT14','RCT1407','RCT1408','RCT183','RCT19','RCT20','RCT29','RCT3','RCT392','RCT4','RCT5','RCT6','RCT8')
for (nname in nnames) mri_lab = mri_lab[!(mri_lab[[nname]] %in% c(-1,-4)),]
mri = mri_lab[, c(colnames(mri), nnames)]

# Data: META+MMSE Features
mri = mri[complete.cases(mri), ]
meta_mri_x = mri
meta_mri_fnames = colnames(mri)[-1]


# write function
writeTGL = function(data, tnames, xname, yname, bname, fnames) {
    for (tname in tnames) {
        y = data[data$VISCODE==tname, yname]
        x = data[data$VISCODE==tname, fnames]
        #print(paste(TD, xname, '_', yname, '_', bname, '_', tname, ':', as.character(dim(x)[1]), sep=''))
        write.csv(y, paste(TD, xname, '_', yname, '_', bname, '_', tname, '_y.csv', sep=''), quote=FALSE, row.names=FALSE);
        write.csv(x, paste(TD, xname, '_', yname, '_', bname, '_', tname, '_x.csv', sep=''), quote=FALSE, row.names=FALSE);
    }
}


# enumerate feature sets
tnames = c('sc', 'm06', 'm12', 'm24', 'm36', 'm48') # 'sc'

mmse = read.csv(paste(FD, 'MMSE.csv', sep=''))
mmse = mmse[mmse$Phase=='ADNI1', ]
adas = read.csv(paste(FD, 'ADASSCORES.csv', sep=''))

for (i in 1:(length(tnames)-1)) {
    mri_tmp = mri_all[mri_all$VISCODE==tnames[i],]
    mri_tmp = mri_tmp[, !colnames(mri_tmp) %in% c('VISCODE')]
    fnames = colnames(mri_tmp)[-1]

    mri_fx = merge(mri_tmp, mri_x)
    meta_mri_fx = merge(mri_tmp, meta_mri_x)
        
    # mmse, start from screening
    mri_mmse = merge(mri_fx, mmse);
    writeTGL(mri_mmse, tnames[(i+1):length(tnames)], 'MRI', 'MMSCORE', tnames[i], c(fnames, mri_fnames))
    meta_mri_mmse = merge(meta_mri_fx, mmse)
    writeTGL(meta_mri_mmse, tnames[(i+1):length(tnames)], 'META_MRI', 'MMSCORE', tnames[i], c(fnames, meta_mri_fnames))

    # adas-cog, start from baseline
    mri_adas = merge(mri_fx, adas)
    writeTGL(mri_adas, tnames[(i+1):length(tnames)], 'MRI', 'TOTAL11', tnames[i], c(fnames, mri_fnames))
    writeTGL(mri_adas, tnames[(i+1):length(tnames)], 'MRI', 'TOTALMOD', tnames[i], c(fnames, mri_fnames))
    meta_mri_adas = merge(meta_mri_fx, adas)
    writeTGL(meta_mri_adas, tnames[(i+1):length(tnames)], 'META_MRI', 'TOTAL11', tnames[i], c(fnames, meta_mri_fnames))
    writeTGL(meta_mri_adas, tnames[(i+1):length(tnames)], 'META_MRI', 'TOTALMOD', tnames[i], c(fnames, meta_mri_fnames))
}
