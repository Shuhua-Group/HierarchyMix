#Author: Shi Zhang
#Date: 06/2021

import pandas as pd
import numpy as np
import math
import argparse

#model1: hierarchical admixture model
#model2: sequential admixture model

#switch number of each haplotype
def get_HapNum(NHap):
	hapNum_dict = {}
	for i in range(NHap):
		haps = Hap_all[i][1]
		n_AB = n_AC = n_AD = n_BC = n_BD = n_CD = 0
		for index in range(0, len(haps)-1):
			if haps.iloc[index,4] == haps.iloc[index+1, 4]:
				if haps.iloc[index,2]==anc_type_dict['A']:
					if haps.iloc[index+1,2]==anc_type_dict['B']:
						n_AB+=1
					elif haps.iloc[index+1,2]==anc_type_dict['C']:
						n_AC+=1
					elif haps.iloc[index+1,2]==anc_type_dict['D']:
						n_AD+=1
				if haps.iloc[index+1,2]==anc_type_dict['A']:
					if haps.iloc[index,2]==anc_type_dict['B']:
						n_AB+=1
					elif haps.iloc[index,2]==anc_type_dict['C']:
						n_AC+=1
					elif haps.iloc[index,2]==anc_type_dict['D']:
						n_AD+=1
				if haps.iloc[index,2]==anc_type_dict['C']:
					if haps.iloc[index+1,2]==anc_type_dict['B']:
						n_BC+=1
					elif haps.iloc[index+1,2]==anc_type_dict['D']:
						n_CD+=1
				if haps.iloc[index+1,2]==anc_type_dict['C']:
					if haps.iloc[index,2]==anc_type_dict['B']:
						n_BC+=1
					elif haps.iloc[index,2]==anc_type_dict['D']:
						n_CD+=1
				if haps.iloc[index,2]==anc_type_dict['B']:
					if haps.iloc[index+1,2]==anc_type_dict['D']:
						n_BD+=1
				if haps.iloc[index,2]==anc_type_dict['D']:
					if haps.iloc[index+1,2]==anc_type_dict['B']:
						n_BD+=1
		hapNum = [n_AB/2, n_AC/2, n_AD/2, n_BC/2, n_BD/2, n_CD/2]
		hapNum_dict[i] = hapNum	
	
	return hapNum_dict

#m, segment Sum and Number of each haplotype
def get_hapmSN(Nhap, cutoff):
	m_dict = {}
	SN_dict = {}
	for i in range(Nhap):
		haps = Hap_all[i][1]
		ASum, BSum, CSum, DSum = 0,0,0,0
		ANum, BNum, CNum, DNum = 0,0,0,0
					
		a1=haps[haps['Anc']==anc_type_dict['A']]
		a2=a1['End']-a1['Start']
		a2 = a2.tolist()
		for k in range(len(a2)):
			if a2[k] > cutoff:
				ASum += a2[k]-cutoff
				ANum += 1

		b1=haps[haps['Anc']==anc_type_dict['B']]
		b2=b1['End']-b1['Start']
		b2 = b2.tolist()
		for k in range(len(b2)):
			if b2[k] > cutoff:
				BSum += b2[k]-cutoff
				BNum += 1

		c1=haps[haps['Anc']==anc_type_dict['C']]
		c2=c1['End']-c1['Start']
		c2 = c2.tolist()
		for k in range(len(c2)):
			if c2[k] > cutoff:
				CSum += c2[k]-cutoff
				CNum += 1

		d1=haps[haps['Anc']==anc_type_dict['D']]		
		d2=d1['End']-d1['Start']
		d2 = d2.tolist()
		for k in range(len(d2)):
			if d2[k] > cutoff:
				DSum += d2[k]-cutoff
				DNum += 1

		sA = ASum+cutoff*ANum
		sB = BSum+cutoff*BNum
		sC = CSum+cutoff*CNum
		sD = DSum+cutoff*DNum
		L = sA+sB+sC+sD
		mA, mB, mC, mD = sA/L, sB/L, sC/L, sD/L

		m_dict[i] = [mA, mB, mC, mD]
		SN_dict[i] = [[ASum, ANum], [BSum, BNum], [CSum, CNum], [DSum, DNum]]

	return m_dict, SN_dict, L

#m of all haplotype in Hap_Index
def getSumM(Hap_Index):
	m = [0,0,0,0]
	for i in Hap_Index:
		for j in range(4):
			m[j] += all_m_dict[i][j]
	sumM = sum(m)
	for i in range(4):
		m[i] /= sumM	
	return m

#u of all haplotype in Hap_Index
def getSumU(Hap_Index, model):
	PindexM = Pindex[model]
	Sum = [0,0,0,0]
	Num = [0,0,0,0]
	for i in Hap_Index:
		for j in range(4):
			Sum[j] += all_SN_dict[i][PindexM[j]][0]
			Num[j] += all_SN_dict[i][PindexM[j]][1]	

	u = [0,0,0,0]
	for i in range(4):
		u[i] = 1/(Sum[i]/Num[i])
	return u

# switch number of all haplotype in Hap_Index
def getSumNKV(Hap_Index):
	NKV = [0,0,0,0,0,0]
	
	for i in Hap_Index:
		for j in range(6):
			NKV[j] += all_hapNum_dict[i][j]
	return NKV

# uKV of model1
def getUKV_M1(mlist, Tlist):
	m1, m2, m3, m4 = mlist[0], mlist[1], mlist[2], mlist[3]
	t1, t2, t3 = Tlist[0], Tlist[1], Tlist[2]

	u = [0,0,0,0,0,0]
	a1 = m1/(m1+m2)
	a2 = 1-a1
	a3 = m3/(m3+m4)
	a4 = m4
	u[0] = m1*a2*(t2-t1)+m1*m2*t1
	u[1] = m1*m3*t1
	u[2] = m1*m4*t1
	u[3] = m2*m3*t1
	u[4] = m2*m4*t1
	u[5] = m3*a4*(t3-t1)+m3*m4*t1
	
	return u

#uKV of model2
def getUKV_M2(mlist, Tlist):
	m1, m2, m3, m4 = mlist[0], mlist[1], mlist[2], mlist[3]
	t1, t2, t3 = Tlist[0], Tlist[1], Tlist[2]

	u = [0,0,0,0,0,0]
	a1 = m1/(m1+m2)
	a2 = 1-a1
	a3 = m3/(1-m4)
	a4 = m4
	u[0] = m1*a2*(t3-t2)+m1*a2*(1-a3)*(t2-t1)+m1*m2*t1
	u[1] = m1*a3*(t2-t1)+m1*m3*t1
	u[2] = m1*m4*t1
	u[3] = m2*a3*(t2-t1)+m2*m3*t1
	u[4] = m2*m4*t1
	u[5] = m3*m4*t1
	
	return u 

# total m under given model
def tranSumM(SumM, model):
	PindexM = Pindex[model]
	m = [0,0,0,0]
	for i in range(4):
		m[i] = SumM[PindexM[i]]
	return m	

# total NKV under given model
def tranSumNKV(SumNKV, model):
	KVIndexM = KVindex[model]
	NKV = [0,0,0,0,0,0]
	for i in range(6):
		NKV[i] = SumNKV[KVIndexM[i]]
	return NKV

# calculate admixture proportion of model1
def M1_admi_pro(mlist):
	m1, m2, m3, m4 = mlist[0], mlist[1], mlist[2], mlist[3]
	a1 = round(m1/(m1+m2),6)
	a2 = round(1-a1,6)
	a3 = round(m3/(m3+m4),6)
	a4 = round(1-a3,6)	
	a5 = round(m1+m2,6)
	a6 = round(1-a5,6)
	return (a1,a2,a3,a4,a5,a6)

# calculate admixture proportion of model2
def M2_admi_pro(mlist):
	m1, m2, m3, m4 = mlist[0], mlist[1], mlist[2], mlist[3]
	a1 = round(m1/(m1+m2),6)
	a2 = round(1-a1, 6)
	a3 = round(m3/(1-m4),6)
	a4 = round(m4,6)
	return (a1,a2,a3,a4)

# estimate admixture time based on the length distribution of ancestral tracts
def getM1_T_len(mlist,ulist):
	m1, m2, m3, m4 = mlist[0], mlist[1], mlist[2], mlist[3]
	u1, u2, u3, u4 = ulist[0], ulist[1], ulist[2], ulist[3]

	t1_AB=(m2*u2-m1*u1)/(m2*(1-m2)-m1*(1-m1))
	t1_CD=(m4*u4-m3*u3)/(m4*(1-m4)-m3*(1-m3))
	t1=(t1_AB+t1_CD)/2
	t2=(m1+m2)*(u1*(1-m2)-u2*(1-m1))/(m2*(1-m2)-m1*(1-m1))+t1
	t3=(m3+m4)*(u3*(1-m4)-u4*(1-m3))/(m4*(1-m4)-m3*(1-m3))+t1

	return (round(t1), round(t2), round(t3))

# estimate admixture time based on the number distribution of ancestral switch points
def getM1_T_trans(m, NKV):
	m1, m2, m3, m4 = m[0], m[1], m[2], m[3]
	a2 = m2/(m1+m2)
	a4 = m4
	t1 = sum(NKV[1:5])/(chrL*Nhap)/((m1+m2)*(m3+m4))	
	t2 = (NKV[0]/(chrL*Nhap)-m1*m2*t1)/(m1*a2)+t1
	t3 = (NKV[5]/(chrL*Nhap)-m3*m4*t1)/(m3*a4)+t1
	return (round(t1), round(t2), round(t3)) 

# estimate admixture time based on the length distribution of ancestral tracts
def getM2_T_len(mlist,ulist):
	m1, m2, m3, m4 = mlist[0], mlist[1], mlist[2], mlist[3]
	u1, u2, u3, u4 = ulist[0], ulist[1], ulist[2], ulist[3]

	t1=u4/(1-m4)
	t2=((1-m4)*u3-(1-m3)*u4)/(1-m3-m4)+t1
	t3_1=(u1*(1-m4)*(1-m3-m4)-u3*(1-m4)*(1-m1-m4)+u4*m4*(m3-m1))/((1-m4)*(1-m1-m3-m4))
	t3_2=(u2*(1-m4)*(1-m3-m4)-u3*(1-m4)*(1-m2-m4)+u4*m4*(m3-m2))/((1-m4)*(1-m2-m3-m4))
	t3=(t3_1+t3_2)/2+t2	

	return (round(t1), round(t2), round(t3))

# estimate admixture time based on the number distribution of ancestral switch points
def getM2_T_trans(m, NKV):
	m1, m2, m3, m4 = m[0], m[1], m[2], m[3]
	a2 = m2/(m1+m2)
	a3 = m3/(1-m4)
	t1 = (NKV[2]+NKV[4]+NKV[5])/(chrL*Nhap)/((m1*m4+m2*m4+m3*m4))
	t2 = ((NKV[1]+NKV[3])/(chrL*Nhap)-(m1*m3+m2*m3)*t1)/(m1*a3+m2*a3)+t1
	t3 = (NKV[0]/(chrL*Nhap)-m1*m2*t1-m1*a2*(1-a3)*(t2-t1))/(m1*a2)+t2
	return (round(t1), round(t2), round(t3))

# calculate likelihood of ancestral switch points
def getLlk(Hap_Index, uKV, model):
	KVindexM = KVindex[model]

	llk = 0
	
	for i in Hap_Index:
		llkInd = 0
		NKV = all_hapNum_dict[i]
		for j in range(6):
			switchN = round(NKV[KVindexM[j]])
			Lambda = uKV[j]*chrL
			llkInd += switchN*math.log(Lambda)-math.log(math.factorial(switchN))-Lambda
		llk += llkInd

	return llk	

def getCI(data_list,k1,k2):
	data_sort= sorted(data_list)
	dataL = data_sort[k1]
	dataR = data_sort[k2]
	return (dataL,dataR)

# select the optimal admixture model
def selectModel(Hap_Index):
	likeli_dict={}
	SumM = getSumM(Hap_Index)
	SumNKV = getSumNKV(Hap_Index)
							
	for model in model1_list:
		SumM_model = tranSumM(SumM, model)	
		SumU_model = getSumU(Hap_Index ,model)	
	
		M1_T_len = getM1_T_len(SumM_model, SumU_model)

		SumNKV_model = tranSumNKV(SumNKV, model)
		M1_T_trans = getM1_T_trans(SumM_model, SumNKV_model)

		if M1_T_trans[2] > M1_T_trans[0] and M1_T_trans[1] > M1_T_trans[0] and M1_T_trans[0] > 0:
			if M1_T_len[2] > M1_T_len[0] and M1_T_len[1] > M1_T_len[0] and M1_T_len[0] > 0:
				uKVlist = getUKV_M1(SumM_model, M1_T_trans)
				locals()['model1_'+model+"_likelihood"] = getLlk(Hap_Index, uKVlist, model)
				likeli_dict["model1_"+model] = locals()["model1_"+model+"_likelihood"]


	for model in model2_list:
		SumM_model = tranSumM(SumM, model)
		SumU_model = getSumU(Hap_Index, model)

		M2_T_len = getM2_T_len(SumM_model, SumU_model)

		SumNKV_model = tranSumNKV(SumNKV, model)
		M2_T_trans = getM2_T_trans(SumM_model, SumNKV_model)
	
		if M2_T_trans[2] > M2_T_trans[1] > M2_T_trans[0] > 0:
			if M2_T_len[2] > M2_T_len[1] > M2_T_len[0] > 0:
				uKVlist = getUKV_M2(SumM_model, M2_T_trans)
				locals()['model2_'+model+'_likelihood']=getLlk(Hap_Index, uKVlist, model)
				likeli_dict['model2_'+model]=locals()['model2_'+model+'_likelihood']

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
	a1_list=[]
	a3_list=[]
	a5_list=[]
	M1N = M2N = 0
	for i in range(nbootstrap):
		Hap_Index_boots = list(np.random.choice(Hap_Index_all, size = Nhap, replace=True))
		boots_SelM = selectModel(Hap_Index_boots)		
		boots_SelM = boots_SelM.split("_")

		m_M = tranSumM(getSumM(Hap_Index_boots), mtype)
		u_M = getSumU(Hap_Index_boots, mtype)

		if mod=='model1':
			if boots_SelM[0]=='model1' and boots_SelM[1]==mtype:
				M1N += 1
					
				a1,a2,a3,a4,a5,a6=M1_admi_pro(m_M)
				t1,t2,t3=getM1_T_len(m_M, u_M)

				a1_list.append(a1)
				a3_list.append(a3)
				a5_list.append(a5)

				t1_list.append(t1)
				t2_list.append(t2)
				t3_list.append(t3)

		elif mod=='model2':
			if boots_SelM[0]=='model2' and boots_SelM[1]==mtype:
				M2N += 1
			
				a1,a2,a3,a4 = M2_admi_pro(m_M)
				t1,t2,t3 = getM2_T_len(m_M, u_M)			
	
				a1_list.append(a1)
				a3_list.append(a3)
		
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
	t1L = sorted(t1_list)[k1]
	t2L = sorted(t2_list)[k1]
	t3L = sorted(t3_list)[k1]
	t1R = sorted(t1_list)[k2]
	t2R = sorted(t2_list)[k2]
	t3R = sorted(t3_list)[k2]
	t_CI={'t1_CI':[t1L,t1R],'t2_CI':[t2L,t2R],'t3_CI':[t3L,t3R]}
	if mod=='model1':
		a1L,a1R=getCI(a1_list,k1,k2)
		a3L,a3R=getCI(a3_list,k1,k2)
		a5L,a5R=getCI(a5_list,k1,k2)
		a_CI={'a1_CI':[a1L,a1R],'a2_CI':[round(1-a1R,6),round(1-a1L,6)],'a3_CI':[a3L,a3R],'a4_CI':[round(1-a3R,6),round(1-a3L,6)],'a5_CI':[a5L,a5R],'a6_CI':[round(1-a5R,6),round(1-a5L,6)]}
		return (t_CI,a_CI,M1N)
	elif mod=='model2':
		a1L,a1R=getCI(a1_list,k1,k2)
		a3L,a3R=getCI(a3_list,k1,k2)
		a_CI={'a1_CI':[a1L,a1R],'a2_CI':[round(1-a1R,6),round(1-a1L,6)],'a3_CI':[a3L,a3R],'a4_CI':[round(1-a3R,6),round(1-a3L,6)]}
		return (t_CI,a_CI,M2N)

def outA():
	m_M = tranSumM(getSumM(Hap_Index_all), all_SelM[1])
	u_M = getSumU(Hap_Index_all, all_SelM[1])

	if all_SelM[0]=='model1':
		a1,a2,a3,a4,a5,a6 = M1_admi_pro(m_M)
		M1_t1,M1_t2,M1_t3=getM1_T_len(m_M, u_M)
			
		fout.write('Best Model:\tHierarchical Admixture Model\n')
		fout.write("\t".join(['', pop1, str(M1_t2)+'(G)', str(a1)])+"\n")
		fout.write("\t".join(['', pop2, str(M1_t2)+'(G)', str(a2)])+"\n")
		fout.write("\t".join(['', pop3, str(M1_t3)+'(G)', str(a3)])+"\n")
		fout.write("\t".join(['', pop4, str(M1_t3)+'(G)', str(a4)])+"\n")
		fout.write("\t".join(['', pop1+'_'+pop2, str(M1_t1)+'(G)', str(a5)])+"\n")
		fout.write("\t".join(['', pop3+'_'+pop4, str(M1_t1)+'(G)', str(a6)])+"\n")

		if Nboots > 0:
			t_CI,a_CI,M1N=bootstrapping(Nboots,ci,cutoff,'model1',all_SelM[1])
			fout.write('--------------------------------------------\n')
			fout.write('Bootstrapping details\n')
			fout.write('Bootstrapping support ratio:'+str(M1N/Nboots*100)+'% ('+str(M1N)+"/"+str(Nboots)+")\n")
			fout.write("\t".join(['', pop1, str(t_CI['t2_CI'])+'(G)', str(a_CI['a1_CI'])])+"\n")
			fout.write("\t".join(['', pop2, str(t_CI['t2_CI'])+'(G)', str(a_CI['a2_CI'])])+"\n")
			fout.write("\t".join(['', pop3, str(t_CI['t3_CI'])+'(G)', str(a_CI['a3_CI'])])+"\n")
			fout.write("\t".join(['', pop4, str(t_CI['t3_CI'])+'(G)', str(a_CI['a4_CI'])])+"\n")
			fout.write("\t".join(['', pop1+'_'+pop2, str(t_CI['t1_CI'])+'(G)', str(a_CI['a5_CI'])])+"\n")
			fout.write("\t".join(['', pop3+'_'+pop4, str(t_CI['t1_CI'])+'(G)', str(a_CI['a6_CI'])])+"\n")
		
		return [a1,a2,a3,a4,a5,a6], [M1_t1,M1_t2,M1_t3]

	else:
		a1,a2,a3,a4 = M2_admi_pro(m_M)
		M2_t1,M2_t2,M2_t3=getM2_T_len(m_M, u_M)

		fout.write("Best Model:\tSequential Admixture Model\n")
		fout.write("\t".join(['', pop1, str(M2_t3)+'(G)', str(a1)])+"\n")
		fout.write("\t".join(['', pop2, str(M2_t3)+'(G)', str(a2)])+"\n")
		fout.write("\t".join(['', pop3, str(M2_t2)+'(G)', str(a3)])+"\n")
		fout.write("\t".join(['', pop4, str(M2_t1)+'(G)', str(a4)])+"\n")

		if Nboots > 0:
			t_CI,a_CI,M2N=bootstrapping(Nboots,ci,cutoff,'model2',all_SelM[1])
			fout.write('--------------------------------------------\n')
			fout.write('Bootstrapping details\n')
			fout.write('Bootstrapping support ratio:'+str(M2N/Nboots*100)+'% ('+str(M2N)+"/"+str(Nboots)+")\n")
			fout.write("\t".join(['', pop1, str(t_CI['t3_CI'])+'(G)', str(a_CI['a1_CI'])])+"\n")
			fout.write("\t".join(['', pop2, str(t_CI['t3_CI'])+'(G)', str(a_CI['a2_CI'])])+"\n")
			fout.write("\t".join(['', pop3, str(t_CI['t2_CI'])+'(G)', str(a_CI['a3_CI'])])+"\n")
			fout.write("\t".join(['', pop4, str(t_CI['t1_CI'])+'(G)', str(a_CI['a4_CI'])])+"\n")

		return [a1,a2,a3,a4], [M2_t1,M2_t2,M2_t3]

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
Nboots=args.bootstrap
ci=args.ci
fout = open(args.output+".txt", 'w')

fseg=pd.read_table(args.input, sep='\t',names=['Start', 'End', 'Anc', 'Hap', 'Chrom'])
fseg=fseg.astype({'Start':'float', 'End':'float'})
fseg=fseg.dropna(axis=0,how='any')
Hap_all = list(fseg.groupby('Hap'))
Nhap = len(Hap_all)
Hap_Index_all = list(range(0, Nhap))

anc_type=list(sorted(set(fseg['Anc'])))
Labellist = ['A', 'B', 'C', 'D']
anc_type_dict={}
for i in range(4):
	anc_type_dict[Labellist[i]] = anc_type[i]

model1_list = ['ABCD', 'ACBD', 'ADBC']
model2_list = ['ABCD','ABDC','ACBD','ACDB','ADBC','ADCB','BCAD','BCDA','BDAC','BDCA','CDAB','CDBA']

# index of ancestral switch points
KVindex = {}
KVindex['ABCD'] = [0, 1, 2, 3, 4, 5]
KVindex['ABDC'] = [0, 2, 1, 4, 3, 5]
KVindex['ACBD'] = [1, 0, 2, 3, 5, 4]
KVindex['ACDB'] = [1, 2, 0, 5, 3, 4]
KVindex['ADBC'] = [2, 0, 1, 4, 5, 3]
KVindex['ADCB'] = [2, 1, 0, 5, 4, 3]
KVindex['BCAD'] = [3, 0, 4, 1, 5, 2]
KVindex['BCDA'] = [3, 4, 0, 5, 1, 2]
KVindex['BDAC'] = [4, 0, 3, 2, 5, 1]
KVindex['BDCA'] = [4, 3, 0, 5, 2, 1]
KVindex['CDAB'] = [5, 1, 3, 2, 4, 0]
KVindex['CDBA'] = [5, 3, 1, 4, 2, 0]

# index of ancestral population
Pindex = {}
Pindex['ABCD'] = [0,1,2,3]
Pindex['ABDC'] = [0,1,3,2]
Pindex['ACBD'] = [0,2,1,3]
Pindex['ACDB'] = [0,2,3,1]
Pindex['ADBC'] = [0,3,1,2]
Pindex['ADCB'] = [0,3,2,1]
Pindex['BCAD'] = [1,2,0,3]
Pindex['BCDA'] = [1,2,3,0]
Pindex['BDAC'] = [1,3,0,2]
Pindex['BDCA'] = [1,3,2,0]
Pindex['CDAB'] = [2,3,0,1]
Pindex['CDBA'] = [2,3,1,0]

all_m_dict, all_SN_dict, chrL = get_hapmSN(Nhap, cutoff)
all_hapNum_dict=get_HapNum(Nhap)

all_SelM = selectModel(Hap_Index_all)

if all_SelM == "none":
	fout.write('No suitable model!\n')
else:
	all_SelM = all_SelM.split('_')
	pop1 = anc_type_dict[all_SelM[1][0]]
	pop2 = anc_type_dict[all_SelM[1][1]]
	pop3 = anc_type_dict[all_SelM[1][2]]
	pop4 = anc_type_dict[all_SelM[1][3]]

	outA()
	
fout.close()

