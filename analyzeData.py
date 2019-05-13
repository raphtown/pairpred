# -*- coding: utf-8 -*-
"""
Created on Tue Oct 09 21:40:36 2012
Check different properties of the data
@author: Afsar
"""

pklpath='./PDBPKL/'
from myPDB import *
import glob
fnames=glob.glob(pklpath+'*l_b.pdb.pkl')
ASA=[]
rASA=[]
RD=[]
CN=[]
UN=[]
DN=[]
PHI=[]
PSI=[]
SS=[]
nRD=[]

for (idx,f) in enumerate(fnames):
    L=myPDB.loader(f)
    lcid=[]
    rcid=[]
    for c in L.stx[0]:
        lcid.append(c.id)
    """
    rname=L.name.split('_')
    rname=rname[0]+'_r_'+rname[2]
    rname=pklpath+rname+'.pdb.pkl'
    
    R=myPDB.loader(rname)
    for c in R.stx[0]:
        rcid.append(c.id)
    cc=list(set(lcid)&set(rcid))
    if len(cc)>0:
        pdb.set_trace()
    
    #L.getStride()
    """
"""
    rASA.append(L.rASA)
    ASA.append(L.ASA)
    RD.append(L.RD[0])
    idxx=idx=np.where(~(np.isnan(L.RD[0])))
    nRD.append(L.RD[0]/L.RD[0][idx].max())
    CN.append(L.CN)
    UN.append(L.UN)
    DN.append(L.DN)
    PHI.append(L.Phi)
    PSI.append(L.Psi)    
    Z=np.nonzero(L.SS>0)
    X=np.empty(L.SS.shape[1]); X.fill(np.nan)
    X[Z[1]]=Z[0]
    SS.append(X)
ZrASA=np.concatenate(rASA)
idx=np.where(~(np.isnan(ZrASA)))
ZrASA=ZrASA[idx]
ZASA=np.concatenate(ASA)[idx]

ZRD=np.concatenate(RD)[idx]
ZnRD=np.concatenate(nRD)[idx]

ZCN=np.concatenate(CN)[idx]
ZUN=np.concatenate(UN)[idx]
ZDN=np.concatenate(DN)[idx]
ZPHI=np.concatenate(PHI)[idx]
ZPSI=np.concatenate(PSI)[idx]
ZSS=np.concatenate(SS)[idx]
#np.where(np.isnan(rASA))
cc=('r.','g.','b.','c.','m.','y.','k.','mx')
for k in range(6):
    idx=ZSS==k
    plt.plot(ZPHI[idx],ZPSI[idx],cc[k])
plt.show()
"""