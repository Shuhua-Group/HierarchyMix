#Author: Shi Zhang
#Date: 06/2021

import pandas as pd
import numpy as np
import math
import argparse

#model1: hierarchical admixture model
#model2: sequential admixture model

#switch number
def num_tran(Hap_Index):
	n_AB=0
	n_AC=n_BC=n_DC=0
	n_AD=0
	n_BA=n_BD=0
	n_CA=n_CB=n_CD=0
	n_DA=n_DB=0

	for i in Hap_Index:
		haps = Hap_all[i][1]
		ifirst = haps.index[0]
		for index in np.linspace(0, len(haps)-2, len(haps)-1):
			index = index+ifirst
			if haps.loc[index,'Anc']==anc_type_dict['A']:
				if haps.loc[index+1,'Anc']==anc_type_dict['B']:
					n_AB+=1
				elif haps.loc[index+1,'Anc']==anc_type_dict['C']:
					n_AC+=1
				elif haps.loc[index+1,'Anc']==anc_type_dict['D']:
					n_AD+=1
			if haps.loc[index+1,'Anc']==anc_type_dict['A']:
				if haps.loc[index,'Anc']==anc_type_dict['B']:
					n_BA+=1
				elif haps.loc[index,'Anc']==anc_type_dict['C']:
					n_CA+=1
				elif haps.loc[index,'Anc']==anc_type_dict['D']:
					n_DA+=1
			if haps.loc[index,'Anc']==anc_type_dict['C']:
				if haps.loc[index+1,'Anc']==anc_type_dict['B']:
					n_CB+=1
				elif haps.loc[index+1,'Anc']==anc_type_dict['D']:
					n_CD+=1
			if haps.loc[index+1,'Anc']==anc_type_dict['C']:
				if haps.loc[index,'Anc']==anc_type_dict['B']:
					n_BC+=1
				elif haps.loc[index,'Anc']==anc_type_dict['D']:
					n_DC+=1
			if haps.loc[index,'Anc']==anc_type_dict['B']:
				if haps.loc[index+1,'Anc']==anc_type_dict['D']:
					n_BD+=1
			if haps.loc[index,'Anc']==anc_type_dict['D']:
				if haps.loc[index+1,'Anc']==anc_type_dict['B']:
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
	Lab=Lba=round(uab*n_len)
	Lac=Lca=round(uac*n_len)
	Lad=Lda=round(uad*n_len)
	Lbc=Lcb=round(ubc*n_len)
	Lbd=Ldb=round(ubd*n_len)
	Lcd=Ldc=round(ucd*n_len)
	return ([Lab,Lac,Lad,Lba,Lbc,Lbd,Lca,Lcb,Lcd,Lda,Ldb,Ldc])

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
	Lab=Lba=round(uab*n_len)
	Lac=Lca=round(uac*n_len)
	Lad=Lda=round(uad*n_len)
	Lbc=Lcb=round(ubc*n_len)
	Lbd=Ldb=round(ubd*n_len)
	Lcd=Ldc=round(ucd*n_len)

	return ([Lab,Lac,Lad,Lba,Lbc,Lbd,Lca,Lcb,Lcd,Lda,Ldb,Ldc])

def model1_T(mA,mB,mC,mD,uA,uB,uC,uD):
	t1_AB=(mB*uB-mA*uA)/(mB*(1-mB)-mA*(1-mA))
	t1_CD=(mD*uD-mC*uC)/(mD*(1-mD)-mC*(1-mC))
	t1=(t1_AB+t1_CD)/2
	t2=(mA+mB)*(uA*(1-mB)-uB*(1-mA))/(mB*(1-mB)-mA*(1-mA))+t1
	t3=(mC+mD)*(uC*(1-mD)-uD*(1-mC))/(mD*(1-mD)-mC*(1-mC))+t1
	return(t1,t2,t3)

def model2_T(mA,mB,mC,mD,uA,uB,uC,uD):
	t1=uD/(1-mD)
	t2=((1-mD)*uC-(1-mC)*uD)/(1-mC-mD)+t1
	t3_1=(uA*(1-mD)*(1-mC-mD)-uC*(1-mD)*(1-mA-mD)+uD*mD*(mC-mA))/((1-mD)*(1-mA-mC-mD))
	t3_2=(uB*(1-mD)*(1-mC-mD)-uC*(1-mD)*(1-mB-mD)+uD*mD*(mC-mB))/((1-mD)*(1-mB-mC-mD))
	if t3_1>=0 and t3_2>=0:
		t3=(t3_1+t3_2)/2+t2
	elif t3_1>0 and t3_2<0:
		t3=t3_1+t2
	elif t3_1<0 and t3_2>0:
		t3=t3_2+t2
	else:
		return None
	if t3>t2>t1:
		return (t1,t2,t3)
	else:
		return None

def getCI(data_list,k1,k2):
	data_sort= sorted(data_list)
	data_lower = round(data_sort[k1],6)
	data_higher = round(data_sort[k2],6)
	return (data_lower,data_higher)

def getm_u(Hap_Index, cutoff, list1, list2):
	Alen = []
	Blen = []
	Clen = []
	Dlen = []

	for i in Hap_Index:
		haps = Hap_all[i][1]
		a1=haps[haps['Anc']==list1[0]]
		a2=a1['End']-a1['Start']
		a2 = a2.tolist()
		for k in range(len(a2)):
			if a2[k] > cutoff:
				Alen.append(a2[k] - cutoff)

		b1=haps[haps['Anc']==list1[1]]
		b2=b1['End']-b1['Start']
		b2 = b2.tolist()
		for k in range(len(b2)):
			if b2[k] > cutoff:
				Blen.append(b2[k] - cutoff)

		c1=haps[haps['Anc']==list1[2]]
		c2=c1['End']-c1['Start']
		c2 = c2.tolist()
		for k in range(len(c2)):
			if c2[k] > cutoff:
				Clen.append(c2[k] - cutoff)
	
		d1=haps[haps['Anc']==list1[3]]
		d2=d1['End']-d1['Start']
		d2 = d2.tolist()
		for k in range(len(d2)):
			if d2[k] > cutoff:
				Dlen.append(d2[k] - cutoff)

	lA = np.mean(Alen)
	uA = 1 / lA

	lB = np.mean(Blen)
	uB = 1 / lB

	lC = np.mean(Clen)
	uC = 1 / lC

	lD = np.mean(Dlen)
	uD = 1 / lD

	sA = (lA+cutoff)*len(Alen)
	sB = (lB+cutoff)*len(Blen)
	sC = (lC+cutoff)*len(Clen)
	sD = (lD+cutoff)*len(Dlen)
	L=sA+sB+sC+sD
	mA, mB, mC, mD = sA/L, sB/L, sC/L, sD/L

	m_dict={}
	m_dict[list2[0]]=mA
	m_dict[list2[1]]=mB
	m_dict[list2[2]]=mC
	m_dict[list2[3]]=mD
	u_dict={}
	u_dict[list2[0]]=uA
	u_dict[list2[1]]=uB
	u_dict[list2[2]]=uC
	u_dict[list2[3]]=uD
	return m_dict, u_dict, L

def selectModel(num_dict, m_dict, u_dict, L):
	likeli_dict={}
	
	for models in model1_list:
		Anc=list(models)
		model1_t_list=model1_T(m_dict[Anc[0]], m_dict[Anc[1]], m_dict[Anc[2]], m_dict[Anc[3]], u_dict[Anc[0]], u_dict[Anc[1]], u_dict[Anc[2]], u_dict[Anc[3]])
		mt11=model1_t_list[0]
		mt211=model1_t_list[1]-model1_t_list[0]
		mt311=model1_t_list[2]-model1_t_list[0]

		if mt311 > 0 and mt211 > 0 and mt11 > 0:
			locals()['model1_lambda_'+models]=model1_lambda(m_dict[Anc[0]],m_dict[Anc[1]],m_dict[Anc[2]],m_dict[Anc[3]],L,mt11,mt211,mt311)
			locals()['model1_'+models+'_likelihood']=getLlk(num_dict[models],locals()['model1_lambda_'+models])
			likeli_dict['model1_'+models]=locals()['model1_'+models+'_likelihood']

	for models in model2_list:
		Anc=list(models)
		model2_t_list = model2_T(m_dict[Anc[0]], m_dict[Anc[1]], m_dict[Anc[2]], m_dict[Anc[3]], u_dict[Anc[0]], u_dict[Anc[1]], u_dict[Anc[2]], u_dict[Anc[3]])
		if model2_t_list != None:
			mt1=model2_t_list[0]
			mt2=model2_t_list[1]
			mt3=model2_t_list[2]

			if mt3 > mt2 > mt1 > 0:
				locals()['model2_lambda_'+models]=model2_lambda(m_dict[Anc[0]],m_dict[Anc[1]],m_dict[Anc[2]],m_dict[Anc[3]],L,mt1,mt2,mt3)
				locals()['model2_'+models+'_likelihood']=getLlk(num_dict[models],locals()['model2_lambda_'+models])
				likeli_dict['model2_'+models]=locals()['model2_'+models+'_likelihood']

	if likeli_dict=={}:
		selection_model = "none"
	else:
		result = sorted(likeli_dict.items(), key=lambda x: x[1], reverse=True)
		selection_model=result[0][0]

	return selection_model	
	
def bootstrapping(nbootstrap, alpha, cutoff, mod, mtype):
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
	for i in range(nbootstrap):
		print("bootstrap "+ str(i)) #
		Hap_Index_boots = list(np.random.choice(Hap_Index_all, size = len(Hap_Index_all), replace=True))

		bootlist1 = [anc_type_dict[anc1], anc_type_dict[anc2], anc_type_dict[anc3], anc_type_dict[anc4]]
		bootlist2 = [anc1, anc2, anc3, anc4]
		boots_m_dict, boots_u_dict, boots_L = getm_u(Hap_Index_boots, cutoff, bootlist1, bootlist2)
		boots_num_dict=num_tran(Hap_Index_boots)

		boots_selection_model = selectModel(boots_num_dict, boots_m_dict, boots_u_dict, boots_L)	
		boots_selection_model = boots_selection_model.split("_")

		if mod=='model1':
			ty = list(mtype)
			if boots_selection_model[0]=='model1' and boots_selection_model[1]==mtype:
				judge_model1+=1
				t1,t2,t3=model1_T(boots_m_dict[ty[0]],boots_m_dict[ty[1]],boots_m_dict[ty[2]],boots_m_dict[ty[3]],boots_u_dict[ty[0]],boots_u_dict[ty[1]],boots_u_dict[ty[2]],boots_u_dict[ty[3]])
				m1,m2,m3,m4,m5,m6=model1_admi_pro(boots_m_dict[ty[0]],boots_m_dict[ty[1]],boots_m_dict[ty[2]],boots_m_dict[ty[3]])

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
			ty=list(mtype)
			if boots_selection_model[0]=='model2' and boots_selection_model[1]==mtype:
				judge_model2+=1
				m1,m2,m3,m4=model2_admi_pro(boots_m_dict[ty[0]],boots_m_dict[ty[1]],boots_m_dict[ty[2]],boots_m_dict[ty[3]])
				t1,t2,t3=model2_T(boots_m_dict[ty[0]],boots_m_dict[ty[1]],boots_m_dict[ty[2]],boots_m_dict[ty[3]],boots_u_dict[ty[0]],boots_u_dict[ty[1]],boots_u_dict[ty[2]],boots_u_dict[ty[3]])

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

def outA(fout):
	if all_selection_model[0]=='model1':
		model1_t1,model1_t2,model1_t3=model1_T(all_m_dict[anc1],all_m_dict[anc2],all_m_dict[anc3],all_m_dict[anc4],all_u_dict[anc1],all_u_dict[anc2],all_u_dict[anc3],all_u_dict[anc4])
		model1_t1=round(model1_t1)
		model1_t2=round(model1_t2)
		model1_t3=round(model1_t3)

		adm1=round(all_m_dict[anc1]/(all_m_dict[anc1]+all_m_dict[anc2]),6)
		adm2=round(1-adm1,6)
		adm3=round(all_m_dict[anc3]/(all_m_dict[anc3]+all_m_dict[anc4]),6)
		adm4=round(1-adm3,6)
		adm5=round(all_m_dict[anc1]+all_m_dict[anc2],6)
		adm6=round(1-adm5,6)

		fout.write('Best Model:\tHierarchical Admixture Model\n')
		fout.write("\t".join(['', anc_type_dict[anc1], str(model1_t2)+'(G)', str(adm1)])+"\n")
		fout.write("\t".join(['', anc_type_dict[anc2], str(model1_t2)+'(G)', str(adm2)])+"\n")
		fout.write("\t".join(['', anc_type_dict[anc3], str(model1_t3)+'(G)', str(adm3)])+"\n")
		fout.write("\t".join(['', anc_type_dict[anc4], str(model1_t3)+'(G)', str(adm4)])+"\n")
		fout.write("\t".join(['', anc_type_dict[anc1]+'_'+anc_type_dict[anc2], str(model1_t1)+'(G)', str(adm5)])+"\n")
		fout.write("\t".join(['', anc_type_dict[anc3]+'_'+anc_type_dict[anc4], str(model1_t1)+'(G)', str(adm6)])+"\n")

		if boots_repeats > 0:
			t_CI,m_CI,model1_count=bootstrapping(boots_repeats,ci,cutoff,'model1',all_selection_model[1])
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
		model2_t1,model2_t2,model2_t3=model2_T(all_m_dict[anc1],all_m_dict[anc2],all_m_dict[anc3],all_m_dict[anc4],all_u_dict[anc1],all_u_dict[anc2],all_u_dict[anc3],all_u_dict[anc4])
		model2_t1=round(model2_t1)
		model2_t2=round(model2_t2)
		model2_t3=round(model2_t3)
		adm1=round(all_m_dict[anc1]/(all_m_dict[anc1]+all_m_dict[anc2]),6)
		adm2=round(1-adm1,6)
		adm3=round(all_m_dict[anc3]/(1-all_m_dict[anc4]),6)
		adm4=round(all_m_dict[anc4],6)

		fout.write("Best Model:\tSequential Admixture Model\n")
		fout.write("\t".join(['', anc_type_dict[anc1], str(model2_t3)+'(G)', str(adm1)])+"\n")
		fout.write("\t".join(['', anc_type_dict[anc2], str(model2_t3)+'(G)', str(adm2)])+"\n")
		fout.write("\t".join(['', anc_type_dict[anc3], str(model2_t2)+'(G)', str(adm3)])+"\n")
		fout.write("\t".join(['', anc_type_dict[anc4], str(model2_t1)+'(G)', str(adm4)])+"\n")

		if boots_repeats > 0:
			t_CI,m_CI,model2_count=bootstrapping(boots_repeats,ci,cutoff,'model2',all_selection_model[1])
			fout.write('--------------------------------------------\n')
			fout.write('Bootstrapping details\n')
			fout.write('Bootstrapping support ratio:'+str(model2_count/boots_repeats*100)+'% ('+str(model2_count)+"/"+str(boots_repeats)+")\n")
			fout.write("\t".join(['', anc_type_dict[anc1], str(t_CI['t3_CI'])+'(G)', str(m_CI['m1_CI'])])+"\n")
			fout.write("\t".join(['', anc_type_dict[anc2], str(t_CI['t3_CI'])+'(G)', str(m_CI['m2_CI'])])+"\n")
			fout.write("\t".join(['', anc_type_dict[anc3], str(t_CI['t2_CI'])+'(G)', str(m_CI['m3_CI'])])+"\n")
			fout.write("\t".join(['', anc_type_dict[anc4], str(t_CI['t1_CI'])+'(G)', str(m_CI['m4_CI'])])+"\n")


parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, required=True, \
										help="Input file name")
parser.add_argument("--lower", type=float, required=False,default=0, \
										help="Lower bound to discard short tracts")
parser.add_argument("--bootstrap", type=int, required=False,default=0, \
										help="Number of bootstrapping")
parser.add_argument("--ci", type=float, required=False, default=0.95, \
										help="The confidence level of bootstrapping confidence interval")
parser.add_argument("--output", type=str, required=False, default="output", \
										help="Prefix of output file")

args = parser.parse_args()
cutoff = args.lower
boots_repeats=args.bootstrap
ci=args.ci
fout = open(args.output+".txt", 'w')

model1_list = ['ABCD', 'ACBD', 'ADBC']
model2_list = ['ABCD','ABDC','ACBD','ACDB','ADBC','ADCB','BCAD','BCDA','BDAC','BDCA','CDAB','CDBA']
		
fseg=pd.read_table(args.input, sep='\t',names=['Start','End','Anc', 'Hap'])
fseg=fseg.astype({'Start':'float','End':'float'})
fseg=fseg.dropna(axis=0,how='any')
anc_type=list(sorted(set(fseg['Anc'])))
anc_type_dict={}
Hap_all = list(fseg.groupby('Hap'))
Hap_Index_all = list(range(0, len(Hap_all)))

Labellist = ['A', 'B', 'C', 'D']
for i in range(4):
	anc_type_dict[Labellist[i]] = anc_type[i]

all_m_dict, all_u_dict, all_L = getm_u(Hap_Index_all, cutoff, anc_type, Labellist)
all_num_dict=num_tran(Hap_Index_all)

all_selection_model = selectModel(all_num_dict, all_m_dict, all_u_dict, all_L)

if all_selection_model == "none":
	fout.write('No suitable model!\n')

all_selection_model=all_selection_model.split('_')
anc1, anc2 = all_selection_model[1][0], all_selection_model[1][1]
anc3, anc4 = all_selection_model[1][2], all_selection_model[1][3]

outA(fout)
	
fout.close()


