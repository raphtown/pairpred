# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 14:09:47 2013
Pymol handling script to go with makeDistanceLinks.py
@author: Afsar
"""

from makeDistanceLinks import readLinkFile
from pymol import cmd
import numpy as np
import pdb
dnames={} #Stores the information about links drawn
"""
dnames[(lmol,lchain,lresi),(rmol,rchain,rresi)]=[dname,lat,lsstr,rat,rsstr,width,rnk,pscore,dmin,disp]
(lmol,lchain,lresi): id of a ligand residue
    lmol: molecule name (string)
    lchain: chain id (char)
    lresi: residue id (integer)
(rmol,rchain,rresi): same as above but for the receptor
dname: name of the distance link created as
    'd'+str(rank)+'_'+lmol+'_'+lchain+str(lresi)+'_'+rmolx+rchain+str(rresi)+'_'+"%1.2f" % pscore+'_'+"%1.2f" % dmin
    where rmolx is empty if lmol==rmol, otherwise it is rmol followd by '_'
lat: identifier of the atom of the ligand residue which is closest to the receptor atom (rat)
rat: receptor atom
lsstr: string representation (pymol) of the ligand residue = '/'+lmol+'//'+lchain+'/'+str(lresi)
rsstr: string representation (pymol) of the ligand residue
width: width of the distance linker\
rnk: rank of the prediction
pscore: prediction score
dmin: minimum distance between the two residues
disp: If true then showLinks will disply that link otherwise it will not be displayed and if it is already being displayed
    then it will br removed from display (not from dnames)    
"""
def sortdnames():
    """
    Return the list of keys in dnames sorted in descending order wrt the prediction socre
    """
    K=dnames.keys()    
    V=dnames.values()
    R=np.array([pscore for (dname,lat,lsstr,rat,rsstr,width,rnk,pscore,dmin,disp) in V])
    si=np.argsort(R)[::-1]
    #print R[si]
    K2=[]
    for i in range(len(si)):
        K2.append(K[si[i]])
        #(dname,lat,lsstr,rat,rsstr,width,rnk,pscore,dmin,disp) = dnames[K[si[i]]]    
        #print pscore
    return K2
    
def map2val(lv,mn,mx,M,m):
    """
    Linear mapping of values in L to the range mn,mx
    
    if np.abs(r)<1e-4:
        r=1.0
        m=M-r    
    """
    r=(M-m)
    v=mn+(mx-mn)*((lv-m)/r)
    #print v
    if v>mx:
        v=mx
    if v<mn:
        v=mn
    return v
def reset():
    """
    Remove links and clear dnames
    """
    global dnames    
    cmd.hide('sticks',"all")
    cmd.hide('lines',"all")
    for dname in dnames:
        cmd.delete(dnames[dname][0])    
    dnames={}
    print "Links Cleared."
def getMinDistAtoms(lsstr,rsstr):
    """
    Returns the atoms that have the minimum distance among the two residues
    passed as input and also the distance between them
    """
    dmin=np.inf
    for at1 in cmd.index(lsstr): 
        for at2 in cmd.index(rsstr): 
            x=cmd.distance('tmp', "%s`%d"%at1, "%s`%d"%at2)            
            if x<dmin:
                dmin=x
                lat=at1
                rat=at2
    cmd.delete('tmp')
    return (lat,rat,dmin)     
def rmvLongEqvLinks(dE):
    """
    Function for disabling the display of equivalent links in a complex
    Two links ((lmol1,lchain1,lresi1),(rmol1,rchain1,rresi1)) and 
    ((lmol2,lchain1,lresi2),(rmol2,rchain1,rresi2)) are considered to be 
    equivalent if 
    (lmol1,lchain1,lresi1) is the same as either of (lmol2,lchain1,lresi2) or (rmol2,rchain1,rresi1)
    OR
    (lmol1,lchain1,lresi1) is the same as either of (lmol2,lchain1,lresi2) or (rmol2,rchain1,rresi1)
    OR
    Both of the above
    If two links are equivalent and they differ significantly in the minimum distance (by more than 1A)
    then only the link with the minimum distance is displayed.
    This is done to clear up the display when predictions for a single chain are to be mapped onto a complex
    in which there are multiple chains corresponding to a single chain in the prediction
    The Equivalence is specified by a dictionary object:
        dE={('3SDL','A'):('3SDL','B'),('3SDL','B'):('3SDL','A'),('3SDL','C'):('3SDL','D'),('3SDL','D'):('3SDL','C')}
    """
    # E is a list of typles of the form ((lmol,lchain),(rmol,rchain))
    # this indicates the equivalence of the two chains
    #print type(dE),dE
    dE={('3SDL','A'):('3SDL','B'),('3SDL','B'):('3SDL','A'),('3SDL','C'):('3SDL','D'),('3SDL','D'):('3SDL','C')}
    for ck in dnames:
        (lmol,lchain,lresi),(rmol,rchain,rresi)=ck
        (_,_,_,_,_,_,_,_,dmin,disp)=dnames[ck]
        if disp:
            kk=[]
            ll=(lmol,lchain) in dE
            rr=(rmol,rchain) in dE
            if ll:                
                (lm,lc)=dE[(lmol,lchain)]
                kk.append(((lm,lc,lresi),(rmol,rchain,rresi)))
                kk.append(((rmol,rchain,rresi),(lm,lc,lresi)))
            if rr:
                (rm,rc)=dE[(rmol,rchain)]
                kk.append(((lmol,lchain,lresi),(rm,rc,rresi)))
                kk.append(((rm,rc,rresi),(lmol,lchain,lresi)))            
            if ll and rr:
                kk.append(((lm,lc,lresi),(rm,rc,rresi)))
                kk.append(((rm,rc,rresi),(lm,lc,lresi)))
            #print kk
            #pdb.set_trace()
            for k in kk:
                if k in dnames:
                    (_,_,_,_,_,_,_,_,kdmin,kdisp)=dnames[k]
                    if kdmin>dmin and (kdmin-dmin)>1.5: 
                        dnames[k][-1]=False                        
                    else:
                        dnames[k][-1]=True 
            
def setDispFlag(dthr=-1):
    """
    sets the display flag of a link based on a distance threshold
    If the minimum distance between the residues forming a link is less 
    than dthr, then the link is displayed, otherwise not
    dthr < 0 implies that all links will be shown
    dthr = 0 implies that no links will be shown
    """    
    dthr=float(dthr)
    if dthr<0:
        dthr=np.inf
    for dname in dnames:
        (dstr,lat,lsstr,rat,rsstr,val,rnk,pscore,dmin,disp)=dnames[dname]
        if dmin<=dthr:
            dnames[dname][-1]=True
        else:
            dnames[dname][-1]=False     
            
def showdnames(dthr,color,mdthr,mcolor,score,mn,mx,mm,MM,dlbl):
    """
    Displays the links that have been selected for display
    sets the color of the links to the specified color
    """
    #calculate the minmum and maximum ranges
    #print dlbl
    if mm == 0 and MM == 0:
        mm=np.inf
        MM=-np.inf
        for dname in dnames:
            (_,_,_,_,_,_,_,pscore,dmin,disp)=dnames[dname]
            if not score:
                stat=-dmin
            else:
                stat=pscore
            if disp:
                if stat<mm:
                    mm=stat
                if stat>MM:
                    MM=stat     
    #plot them
    #print mm,MM
    c=0
    mc=0
    for dname in sortdnames():
        (dstr,lat,lsstr,rat,rsstr,val,rnk,pscore,dmin,disp)=dnames[dname]        
        if disp:
            #print dstr
            if score:
                stat=pscore
            else:
                stat=-dmin
            val=map2val(stat,mn,mx,MM,mm)
            dnames[dname][5]=val 
            #print val
            cmd.distance(dstr, "%s`%d"%lat, "%s`%d"%rat)
            cmd.set('dash_width',val,dstr)
            if dmin<=mdthr:
                cmd.color(mcolor,dstr)
                mc=mc+1
            else:
                cmd.color(color,dstr)
            if not dlbl:
                cmd.hide('labels',dstr)
            else:
                cmd.show('labels',dstr)
            cmd.show('sticks',lsstr)
            cmd.show('sticks',rsstr)
            c=c+1
        else:             
            cmd.delete(dnames[dname][0]) 
    print "Number of Links shown:",c
    print "Number of Links with min. distance less than",mdthr,"A is",mc
def loadLinksList(flist,top=10):
    for f in flist:
        loadLink(f,top)
        
def loadLink(lnkfile,top=10):
    """
    loads a link files and adds the links read to the existing structure
    lnkfile: path to the link file
    top: only the top 'top' links will be read in
    """
    #mNames=cmd.get_names_of_type("object:molecule")    
    top=int(top)
    ((lpdb,lmol),(rpdb,rmol),L)=readLinkFile(lnkfile,top)
    if lpdb==rpdb:
        rmolx=""
    else:
        rmolx=rmol+'_'
    for i in range(len(L)):
        (lchain,lresi)=L[i][0]
        (rchain,rresi)=L[i][1]
        lsstr='/'+lmol+'//'+lchain+'/'+str(lresi)#'c. '+lchain+' and i. '+str(lresi)+' and n. CA'
        rsstr='/'+rmol+'//'+rchain+'/'+str(rresi)#'c. '+rchain+' and i. '+str(rresi)+' and n. CA'
        (lat,rat,dmin)=getMinDistAtoms(lsstr,rsstr)    
        dname='d'+str(i)+'_'+lmol+'_'+lchain+str(lresi)+'_'+rmolx+rchain+str(rresi)+'_'+"%1.2f" % L[i][2]+'_'+"%1.1f" % dmin
        if ((lmol,lchain,lresi),(rmol,rchain,rresi)) not in dnames:
            dnames[(lmol,lchain,lresi),(rmol,rchain,rresi)]=[dname,lat,lsstr,rat,rsstr,1.0,i,L[i][2],dmin,False] #name, rank, prediction score, distance
            
def showLinks(dthr=-1,color='grey',mdthr=6.0,mcolor='red',score=True,mn=1.0,mx=10.0,mm=0,MM=0,dlbl=True,dE={}):
    """
    Display links command
    dthr: distance threshold 
        -1 : show all links irrespective of distance
        otherwise show links such that the min. residue distance is less than dthr
    color: the color for the links
    mdthr: distance threshold for coloring links less than this with mcolor
    mcolor: All residues with minmum distance less than mdthr will be colored in mdthr
    score: if  nonzero then the width of the line increases linearly the prediction score, if 0 it decreases linearly with the minimum distance between the residues    
    mn: minimum line width for the linker
    mx: maximum line width for the linker
    mm: All prediction scores below or equal to mm will be mapped to lines of width mn
    MM: All prediction scores abover or equal to mm will be mapped to lines of width mx
    if both mm and MM are zero, then the displayed link with the minimum score will be 
    displayed with the width of mn and the link with the maximum score will be displayed
    with the width of mx
    dlbl: Whether to show distance values on the linkers (1) or not (0)
    dE: Equivalence dictionary #hard coded for now in rmvLongEqvLinks() 
    """    
    dthr=float(dthr)
    mdthr=float(mdthr)
    mn=float(mn)
    mx=float(mx)
    mm=float(mm)
    MM=float(MM)
    dlbl=int(dlbl)
    score=int(score)
    setDispFlag(dthr)
    if type(dE)==dict:
        rmvLongEqvLinks(dE) 
    showdnames(dthr,color,mdthr,mcolor,score,mn,mx,mm,MM,dlbl)
    
cmd.extend( "loadLink", loadLink)
#cmd.extend( "loadLinksList", loadLink)
cmd.extend( "showLinks", showLinks)
cmd.extend("clrLinks",reset)