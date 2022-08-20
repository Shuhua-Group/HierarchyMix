import pandas as pd
import numpy as np
import math
import argparse

#model1: hierarchical admixture model
#model2: sequential admixture model

#switch number
def num_tran(fseg):
    fseg.index=np.linspace(1,len(fseg),len(fseg))
    n_AB=0
    n_AC=n_BC=n_DC=0
    n_AD=0
    n_BA=n_BD=0
    n_CA=n_CB=n_CD=0
    n_DA=n_DB=0
    
    for index in np.linspace(1,len(fseg)-1,len(fseg)-1):
        if fseg.loc[index,'type']==anc_type_dict['A']:
            if fseg.loc[index+1,'type']==anc_type_dict['B']:
                n_AB+=1
            elif fseg.loc[index+1,'type']==anc_type_dict['C']:
                n_AC+=1
            elif fseg.loc[index+1,'type']==anc_type_dict['D']:
                n_AD+=1
        if fseg.loc[index+1,'type']==anc_type_dict['A']:
            if fseg.loc[index,'type']==anc_type_dict['B']:
                n_BA+=1
            elif fseg.loc[index,'type']==anc_type_dict['C']:
                n_CA+=1
            elif fseg.loc[index,'type']==anc_type_dict['D']:
                n_DA+=1
        if fseg.loc[index,'type']==anc_type_dict['C']:
            if fseg.loc[index+1,'type']==anc_type_dict['B']:
                n_CB+=1
            elif fseg.loc[index+1,'type']==anc_type_dict['D']:
                n_CD+=1
        if fseg.loc[index+1,'type']==anc_type_dict['C']:
            if fseg.loc[index,'type']==anc_type_dict['B']:
                n_BC+=1
            elif fseg.loc[index,'type']==anc_type_dict['D']:
                n_DC+=1
        if fseg.loc[index,'type']==anc_type_dict['B']:
            if fseg.loc[index+1,'type']==anc_type_dict['D']:
                n_BD+=1
        if fseg.loc[index,'type']==anc_type_dict['D']:
            if fseg.loc[index+1,'type']==anc_type_dict['B']:
                n_DB+=1
    num_dict={}
    num_dict['ABCD']=[n_AB,n_AC,n_AD,n_BA,n_BC,n_BD,n_CA,n_CB,n_CD,n_DA,n_DB,n_DC]
    num_dict['ABDC']=[n_AB,n_AD,n_AC,n_BA,n_BD,n_BC,n_DA,n_DB,n_DC,n_CA,n_CB,n_CD]
    num_dict['ACBD']=[n_AC,n_AB,n_AD,n_CA,n_CB,n_CD,n_BA,n_BC,n_BD,n_DA,n_DC,n_DB]
    num_dict['ACDB']=[n_AC,n_AD,n_AB,n_CA,n_CD,n_CB,n_DA,n_DC,n_DB,n_BA,n_BC,n_BD]
    num_dict['ADBC']=[n_AD,n_AB,n_AC,n_DA,n_DB,n_DC,n_BA,n_BD,n_BC,n_CA,n_CD,n_CB]
    num_dict['ADCB']=[n_AD,n_AC,n_AB,n_DA,n_DC,n_DB,n_CA,n_CD,n_CB,n_BA,n_BD,n_BC]
    num_dict['BCAD']=[n_BC,n_BA,n_BD,n_CB,n_CA,n_CD,n_AB,n_AC,n_AD,n_DB,n_DC,n_DA]
    num_dict['BCDA']=[n_BC,n_BD,n_BA,n_CB,n_CD,n_CA,n_DB,n_DC,n_DA,n_AB,n_AC,n_AD]
    num_dict['BDAC']=[n_BD,n_BA,n_BC,n_DB,n_DA,n_DC,n_AB,n_AD,n_AC,n_CB,n_CD,n_CA]
    num_dict['BDCA']=[n_BD,n_BC,n_BA,n_DB,n_DC,n_DA,n_CB,n_CD,n_CA,n_AB,n_AD,n_AC]
    num_dict['CDAB']=[n_CD,n_CA,n_CB,n_DC,n_DA,n_DB,n_AC,n_AD,n_AB,n_BC,n_BD,n_BA]
    num_dict['CDBA']=[n_CD,n_CB,n_CA,n_DC,n_DB,n_DA,n_BC,n_BD,n_BA,n_AC,n_AD,n_AB]
    return num_dict

def getLlk(n_tran, lambda_l):
    llk1 = 0
    llk2 = 0
    for i in range(len(n_tran)):
        llk1+=n_tran[i]*math.log(lambda_l[i])
        llk2+=math.log(math.factorial(n_tran[i]))
    llk=llk1-llk2-sum(lambda_l)
    return llk

def model1_admi_pro(a,b,c,d):
    aA=a/(a+b)
    aB=b/(a+b)
    aC=c/(c+d)
    aD=d/(c+d)
    aE=a+b
    aF=1-(a+b)

    return [aA,aB,aC,aD,aE,aF]

def model1_lambda(mA,mB,mC,mD,n_len,t1,t2_AB,t3_CD):
    aA=mA/(mA+mB)
    aB=1-aA
    aC=mC/(mC+mD)
    aD=1-aC
    uab=mA*aB*t2_AB+mA*mB*t1
    uac=mA*mC*t1
    uad=mA*mD*t1
    ubc=mB*mC*t1
    ubd=mB*mD*t1
    ucd=mC*aD*t3_CD+mC*mD*t1
    lambda_ab=lambda_ba=round(uab*n_len)
    lambda_ac=lambda_ca=round(uac*n_len)
    lambda_ad=lambda_da=round(uad*n_len)
    lambda_bc=lambda_cb=round(ubc*n_len)
    lambda_bd=lambda_db=round(ubd*n_len)
    lambda_cd=lambda_dc=round(ucd*n_len)
    
    return ([lambda_ab,lambda_ac,lambda_ad,lambda_ba,lambda_bc,lambda_bd,lambda_ca,lambda_cb,lambda_cd,lambda_da,lambda_db,lambda_dc])

def model2_admi_pro(a,b,c,d):
    aA=a/(a+b)
    aB=b/(a+b)
    aC=c/(1-d)
    aD=d
    return [aA,aB,aC,aD]

def model2_lambda(mA,mB,mC,mD,n_len,t1,t2,t3):
    aA=mA/(mA+mB)
    aB=1-aA
    aC=mC/(1-mD)
    aD=mD
    uab=mA*aB*(t3-t2)+mA*aB*(1-aC)*(t2-t1)+mA*mB*t1
    uac=mA*aC*(t2-t1)+mA*mC*t1
    uad=mA*mD*t1
    ubc=mB*aC*(t2-t1)+mB*mC*t1
    ubd=mB*mD*t1
    ucd=mC*mD*t1
    lambda_ab=lambda_ba=round(uab*n_len)
    lambda_ac=lambda_ca=round(uac*n_len)
    lambda_ad=lambda_da=round(uad*n_len)
    lambda_bc=lambda_cb=round(ubc*n_len)
    lambda_bd=lambda_db=round(ubd*n_len)
    lambda_cd=lambda_dc=round(ucd*n_len)    
    
    return ([lambda_ab,lambda_ac,lambda_ad,lambda_ba,lambda_bc,lambda_bd,lambda_ca,lambda_cb,lambda_cd,lambda_da,lambda_db,lambda_dc])

def model1_T(mA,mB,mC,mD,uA,uB,uC,uD):
    model1_t1_AB=(mB*uB-mA*uA)/(mB*(1-mB)-mA*(1-mA))
    model1_t1_CD=(mD*uD-mC*uC)/(mD*(1-mD)-mC*(1-mC))
    model1_t1=(model1_t1_AB+model1_t1_CD)/2
    model1_t2=(mA+mB)*(uA*(1-mB)-uB*(1-mA))/(mB*(1-mB)-mA*(1-mA))+model1_t1
    model1_t3=(mC+mD)*(uC*(1-mD)-uD*(1-mC))/(mD*(1-mD)-mC*(1-mC))+model1_t1
    return(model1_t1,model1_t2,model1_t3)

def model2_T(mA,mB,mC,mD,uA,uB,uC,uD):
    model2_t1=uD/(1-mD)
    model2_t2=((1-mD)*uC-(1-mC)*uD)/(1-mC-mD)+model2_t1
    model2_t3_1=(uA*(1-mD)*(1-mC-mD)-uC*(1-mD)*(1-mA-mD)+uD*mD*(mC-mA))/((1-mD)*(1-mA-mC-mD))
    model2_t3_2=(uB*(1-mD)*(1-mC-mD)-uC*(1-mD)*(1-mB-mD)+uD*mD*(mC-mB))/((1-mD)*(1-mB-mC-mD))
    if model2_t3_1>=0 and model2_t3_2>=0:
        model2_t3=(model2_t3_1+model2_t3_2)/2+model2_t2
    elif model2_t3_1>0 and model2_t3_2<0:
        model2_t3=model2_t3_1+model2_t2
    elif model2_t3_1<0 and model2_t3_2>0:
        model2_t3=model2_t3_2+model2_t2
    else:
        return None
    if model2_t3>model2_t2>model2_t1:
        return (model2_t1,model2_t2,model2_t3)
    else:
        return None

def getCI(data_list,k1,k2):
    data_sort= sorted(data_list)
    data_lower = round(data_sort[k1],6)
    data_higher = round(data_sort[k2],6)
    return (data_lower,data_higher)

def bootstrapping(f, nbootstrap, alpha, cutoff, mod,ty):
    t1_list=[]
    t2_list=[]
    t3_list=[]
    m1_list=[]
    m2_list=[]
    m3_list=[]
    m4_list=[]
    m5_list=[]
    m6_list=[]
    judge_model1=0
    judge_model2=0
    judge_none=0
    for i in range(nbootstrap):
        fb=pd.DataFrame()
        fb=f.sample(n=len(f),replace=True)
        a1=fb[fb['type']==anc_type_dict[anc1]]
        a2=a1['End']-a1['Start']
        a2 = a2.tolist()
        a3 = []
        for k in range(len(a2)):
            if a2[k] > cutoff:
                a3.append(a2[k] - cutoff)
        lA = np.mean(a3)
        uA = 1 / lA
        
        b1=fb[fb['type']==anc_type_dict[anc2]]
        b2=b1['End']-b1['Start']
        b2 = b2.tolist()
        b3 = []
        for k in range(len(b2)):
            if b2[k] > cutoff:
                b3.append(b2[k] - cutoff)
        lB = np.mean(b3)
        uB = 1 / lB
        
        c1=fb[fb['type']==anc_type_dict[anc3]]
        c2=c1['End']-c1['Start']
        c2 = c2.tolist()
        c3 = []
        for k in range(len(c2)):
            if c2[k] > cutoff:
                c3.append(c2[k] - cutoff)
        lC = np.mean(c3)
        uC = 1 / lC
        
        d1=fb[fb['type']==anc_type_dict[anc4]]
        d2=d1['End']-d1['Start']
        d2 = d2.tolist()
        d3 = []
        for k in range(len(d2)):
            if d2[k] > cutoff:
                d3.append(d2[k] - cutoff)
        lD = np.mean(d3)
        uD = 1 / lD
        length=sum(a2+b2+c2+d2)
        mA=sum(a2)/sum(a2+b2+c2+d2)
        mB=sum(b2)/sum(a2+b2+c2+d2)
        mC=sum(c2)/sum(a2+b2+c2+d2)
        mD=sum(d2)/sum(a2+b2+c2+d2)
        
        likeli_dict={}
        boots_m_dict={}
        boots_m_dict[anc1]=mA
        boots_m_dict[anc2]=mB
        boots_m_dict[anc3]=mC
        boots_m_dict[anc4]=mD
        boots_u_dict={}
        boots_u_dict[anc1]=uA
        boots_u_dict[anc2]=uB
        boots_u_dict[anc3]=uC
        boots_u_dict[anc4]=uD
        sample_num_dict=num_tran(fb)

        for models in model1_list:
            ances=list(models)
            model1_t_list=model1_T(boots_m_dict[ances[0]], boots_m_dict[ances[1]], boots_m_dict[ances[2]], boots_m_dict[ances[3]], boots_u_dict[ances[0]], boots_u_dict[ances[1]], boots_u_dict[ances[2]], boots_u_dict[ances[3]])
            mt11=model1_t_list[0]
            mt211=model1_t_list[1]-model1_t_list[0]
            mt311=model1_t_list[2]-model1_t_list[0]
            
            if mt311 > 0 and mt211 > 0 and mt11 > 0:
                locals()['model1_lambda_'+models]=model1_lambda(boots_m_dict[ances[0]],boots_m_dict[ances[1]],boots_m_dict[ances[2]],boots_m_dict[ances[3]],length,mt11,mt211,mt311)
                locals()['model1_'+models+'_likelihood']=getLlk(sample_num_dict[models],locals()['model1_lambda_'+models])
                likeli_dict['model1_'+models]=locals()['model1_'+models+'_likelihood']

        for models in model2_list:
            ances=list(models)
            model2_t_list = model2_T(boots_m_dict[ances[0]], boots_m_dict[ances[1]], boots_m_dict[ances[2]], boots_m_dict[ances[3]], boots_u_dict[ances[0]], boots_u_dict[ances[1]], boots_u_dict[ances[2]], boots_u_dict[ances[3]])
            if model2_t_list != None:
                mt1=model2_t_list[0]
                mt2=model2_t_list[1]
                mt3=model2_t_list[2]
                if mt3 > mt2 > mt1 > 0:
                    locals()['model2_lambda_'+models]=model2_lambda(boots_m_dict[ances[0]],boots_m_dict[ances[1]],boots_m_dict[ances[2]],boots_m_dict[ances[3]],length,mt1,mt2,mt3)
                    locals()['model2_'+models+'_likelihood']=getLlk(sample_num_dict[models],locals()['model2_lambda_'+models])
                    likeli_dict['model2_'+models]=locals()['model2_'+models+'_likelihood']

        if likeli_dict=={}:
            judge_none+=1
            continue
        else:
            result = sorted(likelihood_dict.items(), key=lambda x: x[1], reverse=True)
            boots_selection_model=result[0][0].split('_')
        if mod=='model1':
            types=list(ty)
            if boots_selection_model[0]=='model1' and boots_selection_model[1]==ty:
                judge_model1+=1
                t1,t2,t3=model1_T(boots_m_dict[types[0]],boots_m_dict[types[1]],boots_m_dict[types[2]],boots_m_dict[types[3]],boots_u_dict[types[0]],boots_u_dict[types[1]],boots_u_dict[types[2]],boots_u_dict[types[3]])
                m1,m2,m3,m4,m5,m6=model1_admi_pro(boots_m_dict[types[0]],boots_m_dict[types[1]],boots_m_dict[types[2]],boots_m_dict[types[3]])
                m1_list.append(m1)
                m2_list.append(m2)
                m3_list.append(m3)
                m4_list.append(m4)
                m5_list.append(m5)
                m6_list.append(m6)
                
                t1_list.append(t1)
                t2_list.append(t2)
                t3_list.append(t3)
                
        elif mod=='model2':
            types=list(ty)
            if boots_selection_model[0]=='model2' and boots_selection_model[1]==ty:
                judge_model2+=1
                m1,m2,m3,m4=model2_admi_pro(boots_m_dict[types[0]],boots_m_dict[types[1]],boots_m_dict[types[2]],boots_m_dict[types[3]])
                t1,t2,t3=model2_T(boots_m_dict[types[0]],boots_m_dict[types[1]],boots_m_dict[types[2]],boots_m_dict[types[3]],boots_u_dict[types[0]],boots_u_dict[types[1]],boots_u_dict[types[2]],boots_u_dict[types[3]])
                m1_list.append(m1)
                m2_list.append(m2)
                m3_list.append(m3)
                m4_list.append(m4)
                
                t1_list.append(t1)
                t2_list.append(t2)
                t3_list.append(t3)

    a = 1 - alpha
    count=len(t1_list)
    k1 = int(count * a / 2)-1
    if k1<0:
        k1=0
    k2 = int(count * (1 - a / 2))-1
    if k2==count:
        k2=count-1
    t1_sort= sorted(t1_list)
    t2_sort= sorted(t2_list)
    t3_sort= sorted(t3_list)
    t1_lower = round(t1_sort[k1])
    t2_lower = round(t2_sort[k1])
    t3_lower = round(t3_sort[k1])
    t1_higher = round(t1_sort[k2])
    t2_higher = round(t2_sort[k2])
    t3_higher = round(t3_sort[k2])
    t_CI={'t1_CI':[t1_lower,t1_higher],'t2_CI':[t2_lower,t2_higher],'t3_CI':[t3_lower,t3_higher]}
    if mod=='model1':
        m1_lower,m1_higher=getCI(m1_list,k1,k2)
        m2_lower,m2_higher=getCI(m2_list,k1,k2)
        m3_lower,m3_higher=getCI(m3_list,k1,k2)
        m4_lower,m4_higher=getCI(m4_list,k1,k2)
        m5_lower,m5_higher=getCI(m5_list,k1,k2)
        m6_lower,m6_higher=getCI(m6_list,k1,k2)
        m_CI={'m1_CI':[m1_lower,m1_higher],'m2_CI':[m2_lower,m2_higher],'m3_CI':[m3_lower,m3_higher],'m4_CI':[m4_lower,m4_higher],'m5_CI':[m5_lower,m5_higher],'m6_CI':[m6_lower,m6_higher]}
        return (t_CI,m_CI,judge_model1)
    elif mod=='model2':
        m1_lower,m1_higher=getCI(m1_list,k1,k2)
        m2_lower,m2_higher=getCI(m2_list,k1,k2)
        m3_lower,m3_higher=getCI(m3_list,k1,k2)
        m4_lower,m4_higher=getCI(m4_list,k1,k2)
        m_CI={'m1_CI':[m1_lower,m1_higher],'m2_CI':[m2_lower,m2_higher],'m3_CI':[m3_lower,m3_higher],'m4_CI':[m4_lower,m4_higher]}
        return (t_CI,m_CI,judge_model2)


parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, required=True, \
                                        help="Input file name")
parser.add_argument("--lower", type=float, required=False,default=0, \
                                        help="Lower bound to discard short tracts")
parser.add_argument("--bootstrap", type=int, required=False,default=0, \
                                        help="Number of bootstrapping")
parser.add_argument("--ci", type=float, required=False, default=0.95, \
                                        help="The confidence level of bootstrapping	confidence interval")
parser.add_argument("--left", type=str, required=False, default=None, \
                                        help="The specified left part of admixture model, e.g., Anc1_Anc2")
parser.add_argument("--right", type=str, required=False, default=None, \
                                        help="The specified right part of admixture model, e.g., Anc3_Anc4")
parser.add_argument("--output", type=str, required=False, default="output", \
										help="Prefix of output file")

args = parser.parse_args()

cutoff = args.lower
boots_repeats=args.bootstrap
ci=args.ci
fout = open(args.output+".txt", 'w')

model1_t_list=[]
if args.left!=None:
    model1_list=['ABCD','ACBD','ADBC']
else:
    model1_list=['ABCD','ABDC','ACBD','ACDB','ADBC','ADCB','BCAD','BCDA','BDAC','BDCA','CDAB','CDBA']
model2_list=['ABCD','ABDC','ACBD','ACDB','ADBC','ADCB','BCAD','BCDA','BDAC','BDCA','CDAB','CDBA']

likelihood_dict={}
infe_l_dict={}
fseg=pd.read_table(args.input, sep='\t',names=['Start','End','type'])
fseg=fseg.astype({'Start':'float','End':'float'})
fseg=fseg.dropna(axis=0,how='any')
anc_type=list(sorted(set(fseg['type'])))
anc_type_dict={}

if args.left==None:
    a1=fseg[fseg['type']==anc_type[0]]
    a2=a1['End']-a1['Start']
    a2 = a2.tolist()
    a3 = []
    for k in range(len(a2)):
        if a2[k] > cutoff:
            a3.append(a2[k] - cutoff)
    lA = np.mean(a3)
    uA = 1 / lA
    
    b1=fseg[fseg['type']==anc_type[1]]
    b2=b1['End']-b1['Start']
    b2 = b2.tolist()
    b3 = []
    for k in range(len(b2)):
        if b2[k] > cutoff:
            b3.append(b2[k] - cutoff)
    lB = np.mean(b3)
    uB = 1 / lB
    
    c1=fseg[fseg['type']==anc_type[2]]
    c2=c1['End']-c1['Start']
    c2 = c2.tolist()
    c3 = []
    for k in range(len(c2)):
        if c2[k] > cutoff:
            c3.append(c2[k] - cutoff)
    lC = np.mean(c3)
    uC = 1 / lC
    
    d1=fseg[fseg['type']==anc_type[3]]
    d2=d1['End']-d1['Start']
    d2 = d2.tolist()
    d3 = []
    for k in range(len(d2)):
        if d2[k] > cutoff:
            d3.append(d2[k] - cutoff)
    lD = np.mean(d3)
    uD = 1 / lD
    
    mA=sum(a2)/sum(a2+b2+c2+d2)
    mB=sum(b2)/sum(a2+b2+c2+d2)
    mC=sum(c2)/sum(a2+b2+c2+d2)
    mD=sum(d2)/sum(a2+b2+c2+d2)
    n_len=sum(a2+b2+c2+d2)
    m=[mA,mB,mC,mD]
    u=[uA,uB,uC,uD]
    
    m_dict={}
    m_dict['A']=mA
    m_dict['B']=mB
    m_dict['C']=mC
    m_dict['D']=mD
    
    u_dict={}
    u_dict['A']=uA
    u_dict['B']=uB
    u_dict['C']=uC
    u_dict['D']=uD
    
    anc_type_dict['A']=anc_type[0]
    anc_type_dict['B']=anc_type[1]
    anc_type_dict['C']=anc_type[2]
    anc_type_dict['D']=anc_type[3]
else:
    left_anc=args.left.split('_')
    right_anc=args.right.split('_')
    a1=fseg[fseg['type']==left_anc[0]]
    a2=a1['End']-a1['Start']
    a2 = a2.tolist()
    a3 = []
    for k in range(len(a2)):
        if a2[k] > cutoff:
            a3.append(a2[k] - cutoff)
    lA = np.mean(a3)
    uA = 1 / lA
    
    b1=fseg[fseg['type']==left_anc[1]]
    b2=b1['End']-b1['Start']
    b2 = b2.tolist()
    b3 = []
    for k in range(len(b2)):
        if b2[k] > cutoff:
            b3.append(b2[k] - cutoff)
    lB = np.mean(b3)
    uB = 1 / lB
    
    c1=fseg[fseg['type']==right_anc[0]]
    c2=c1['End']-c1['Start']
    c2 = c2.tolist()
    c3 = []
    for k in range(len(c2)):
        if c2[k] > cutoff:
            c3.append(c2[k] - cutoff)
    lC = np.mean(c3)
    uC = 1 / lC
    
    d1=fseg[fseg['type']==right_anc[1]]
    d2=d1['End']-d1['Start']
    d2 = d2.tolist()
    d3 = []
    for k in range(len(d2)):
        if d2[k] > cutoff:
            d3.append(d2[k] - cutoff)
    lD = np.mean(d3)
    uD = 1 / lD
    
    n_len=sum(a2+b2+c2+d2)
    mA=sum(a2)/n_len
    mB=sum(b2)/n_len
    mC=sum(c2)/n_len
    mD=sum(d2)/n_len
    m_dict={}
    u_dict={}

    if mA>mB:
        m_dict['A']=mA
        m_dict['B']=mB
        u_dict['A']=uA
        u_dict['B']=uB
        anc_type_dict['A']=left_anc[0]
        anc_type_dict['B']=left_anc[1]
    else:
        m_dict['A']=mB
        m_dict['B']=mA
        u_dict['A']=uB
        u_dict['B']=uA
        anc_type_dict['A']=left_anc[1]
        anc_type_dict['B']=left_anc[0]
    if mC>mD:
        m_dict['C']=mC
        m_dict['D']=mD
        u_dict['C']=uC
        u_dict['D']=uD
        anc_type_dict['C']=right_anc[0]
        anc_type_dict['D']=right_anc[1]
    else:
        m_dict['C']=mD
        m_dict['D']=mC
        u_dict['C']=uD
        u_dict['D']=uC
        anc_type_dict['C']=right_anc[1]
        anc_type_dict['D']=right_anc[0]
    #print(anc_type_dict)

model_num_dict=num_tran(fseg)

for model in model1_list:
    ancestry=list(model)
    model1_t_list = model1_T(m_dict[ancestry[0]], m_dict[ancestry[1]], m_dict[ancestry[2]], m_dict[ancestry[3]], u_dict[ancestry[0]], u_dict[ancestry[1]], u_dict[ancestry[2]], u_dict[ancestry[3]])
    mt1=model1_t_list[0]
    mt21=model1_t_list[1]-model1_t_list[0]
    mt31=model1_t_list[2]-model1_t_list[0]
    if mt31 > 0 and mt21 > 0 and mt1 > 0:
        locals()['model1_lambda_'+model]=model1_lambda(m_dict[ancestry[0]],m_dict[ancestry[1]],m_dict[ancestry[2]],m_dict[ancestry[3]],n_len,mt1,mt21,mt31)
        locals()['model1_'+model+'_likelihood']=getLlk(model_num_dict[model],locals()['model1_lambda_'+model])
        likelihood_dict['model1_'+model]=locals()['model1_'+model+'_likelihood']

for model in model2_list:
    ancestry=list(model)
    model2_t_list = model2_T(m_dict[ancestry[0]], m_dict[ancestry[1]], m_dict[ancestry[2]], m_dict[ancestry[3]], u_dict[ancestry[0]], u_dict[ancestry[1]], u_dict[ancestry[2]], u_dict[ancestry[3]])
    if model2_t_list != None:
        mt1=model2_t_list[0]
        mt2=model2_t_list[1]
        mt3=model2_t_list[2]
        if mt3 > mt2 > mt1 > 0:
            locals()['model2_lambda_'+model]=model2_lambda(m_dict[ancestry[0]],m_dict[ancestry[1]],m_dict[ancestry[2]],m_dict[ancestry[3]],n_len,mt1,mt2,mt3)
            locals()['model2_'+model+'_likelihood']=getLlk(model_num_dict[model],locals()['model2_lambda_'+model])
            likelihood_dict['model2_'+model]=locals()['model2_'+model+'_likelihood']

result = sorted(likelihood_dict.items(), key=lambda x: x[1], reverse=True)

if result==[]:
    fout.write('No suitable model!\n')

selection_model=result[0][0].split('_')
anc1=selection_model[1][0]
anc2=selection_model[1][1]
anc3=selection_model[1][2]
anc4=selection_model[1][3]
if selection_model[0]=='model1':
    if boots_repeats > 0:
        t_CI,m_CI,model1_count=bootstrapping(fseg,boots_repeats,ci,cutoff,'model1',selection_model[1])
    
    model1_t1,model1_t2,model1_t3=model1_T(m_dict[anc1],m_dict[anc2],m_dict[anc3],m_dict[anc4],u_dict[anc1],u_dict[anc2],u_dict[anc3],u_dict[anc4])
    model1_t1=round(model1_t1)
    model1_t2=round(model1_t2)
    model1_t3=round(model1_t3)

    adm1=round(m_dict[anc1]/(m_dict[anc1]+m_dict[anc2]),6)
    adm2=round(1-adm1,6)
    adm3=round(m_dict[anc3]/(m_dict[anc3]+m_dict[anc4]),6)
    adm4=round(1-adm3,6)
    adm5=round(m_dict[anc1]+m_dict[anc2],6)
    adm6=round(1-adm5,6)
    fout.write('Best Model:\tHierarchical Admixture Model\n')
    fout.write("\t".join(['', anc_type_dict[anc1], str(model1_t2)+'(G)', str(adm1)])+"\n")
    fout.write("\t".join(['', anc_type_dict[anc2], str(model1_t2)+'(G)', str(adm2)])+"\n")
    fout.write("\t".join(['', anc_type_dict[anc3], str(model1_t3)+'(G)', str(adm3)])+"\n")
    fout.write("\t".join(['', anc_type_dict[anc4], str(model1_t3)+'(G)', str(adm4)])+"\n")
    fout.write("\t".join(['', anc_type_dict[anc1]+'_'+anc_type_dict[anc2], str(model1_t1)+'(G)', str(adm5)])+"\n")
    fout.write("\t".join(['', anc_type_dict[anc3]+'_'+anc_type_dict[anc4], str(model1_t1)+'(G)', str(adm6)])+"\n")

    if boots_repeats > 0:
        fout.write('--------------------------------------------\n')
        fout.write('Bootstrapping details\n')
        fout.write('Bootstrapping support ratio:'+str(model1_count/boots_repeats*100)+'% ('+str(model1_count)+"/"+str(boots_repeats)+")\n")
        fout.write("\t".join(['', anc_type_dict[anc1], str(t_CI['t2_CI'])+'(G)', str(m_CI['m1_CI'])])+"\n")
        fout.write("\t".join(['', anc_type_dict[anc2], str(t_CI['t2_CI'])+'(G)', str(m_CI['m2_CI'])])+"\n")
        fout.write("\t".join(['', anc_type_dict[anc3], str(t_CI['t3_CI'])+'(G)', str(m_CI['m3_CI'])])+"\n")
        fout.write("\t".join(['', anc_type_dict[anc4], str(t_CI['t3_CI'])+'(G)', str(m_CI['m4_CI'])])+"\n")
        fout.write("\t".join(['', anc_type_dict[anc1]+'_'+anc_type_dict[anc2], str(t_CI['t1_CI'])+'(G)', str(m_CI['m5_CI'])])+"\n")
        fout.write("\t".join(['', anc_type_dict[anc3]+'_'+anc_type_dict[anc4], str(t_CI['t1_CI'])+'(G)', str(m_CI['m6_CI'])])+"\n")

else:
    if boots_repeats > 0:
        t_CI,m_CI,model2_count=bootstrapping(fseg,boots_repeats,ci,cutoff,'model2',selection_model[1])
    
    model2_t1,model2_t2,model2_t3=model2_T(m_dict[anc1],m_dict[anc2],m_dict[anc3],m_dict[anc4],u_dict[anc1],u_dict[anc2],u_dict[anc3],u_dict[anc4])
    model2_t1=round(model2_t1)
    model2_t2=round(model2_t2)
    model2_t3=round(model2_t3)
    adm1=round(m_dict[anc1]/(m_dict[anc1]+m_dict[anc2]),6)
    adm2=round(1-adm1,6)
    adm3=round(m_dict[anc3]/(1-m_dict[anc4]),6)
    adm4=round(m_dict[anc4],6)
    
    fout.write("Best Model:\tSequential Admixture Model\n")
    fout.write("\t".join(['', anc_type_dict[anc1], str(model2_t3)+'(G)', str(adm1)])+"\n")
    fout.write("\t".join(['', anc_type_dict[anc2], str(model2_t3)+'(G)', str(adm2)])+"\n")
    fout.write("\t".join(['', anc_type_dict[anc3], str(model2_t2)+'(G)', str(adm3)])+"\n")
    fout.write("\t".join(['', anc_type_dict[anc4], str(model2_t1)+'(G)', str(adm4)])+"\n")

    if boots_repeats > 0:
        fout.write('--------------------------------------------\n')
        fout.write('Bootstrapping details\n')
        fout.write('Bootstrapping support ratio:'+str(model2_count/boots_repeats*100)+'% ('+str(model2_count)+"/"+str(boots_repeats)+")\n")
        fout.write("\t".join(['', anc_type_dict[anc1], str(t_CI['t3_CI'])+'(G)', str(m_CI['m1_CI'])])+"\n")	
        fout.write("\t".join(['', anc_type_dict[anc2], str(t_CI['t3_CI'])+'(G)', str(m_CI['m2_CI'])])+"\n")
        fout.write("\t".join(['', anc_type_dict[anc3], str(t_CI['t2_CI'])+'(G)', str(m_CI['m3_CI'])])+"\n")
        fout.write("\t".join(['', anc_type_dict[anc4], str(t_CI['t1_CI'])+'(G)', str(m_CI['m4_CI'])])+"\n")

fout.close()

