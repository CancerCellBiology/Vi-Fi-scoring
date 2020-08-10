# -*- coding: utf-8 -*-
"""
Code developed in the laboratory of Cancer Cell Biology,
Engelhardt Institute of Molecular Biology, Moscow, Russia.
Not for commercial use.
If you have any questions about the code or its use contact:
lebedevtd@gmai.com
vr.elmira@gmail.com 
"""

import pandas as pd
import os

def Score(df, score):
    S_up= df.groupby(['change']).sum().at['up', score]/df[df['change']=='up'].shape[0]
    S_down= df.groupby(['change']).sum().at['down', score]/df[df['change']=='down'].shape[0]
    S= S_up-S_down
    return S

"""define working directory"""
os.chdir('C:\\Vi-Fi scoring')
input_file=input('Drug list:')
df_user= pd.read_excel(input_file+'.xlsx')
drugs=list(dict.fromkeys(df_user['Name'].dropna()))
other_drugs=list(dict.fromkeys(df_user['Name'].dropna()))
df_score=pd.DataFrame()
for drug1 in drugs:
    FS_list=[]
    FS1_list=[]
    FS2_list=[]
    FS_change_list=[]
    VS_list=[]
    VS1_list=[]
    VS2_list=[]
    VS_change_list=[]
    names_list=[]
    other_drugs.remove(drug1)
    for drug2 in other_drugs:        
        dfV1= pd.read_excel('Vi-Fi Signatures\\'+drug1+'_Viral_sig.xlsx')
        dfV2= pd.read_excel('Vi-Fi Signatures\\'+drug2+'_Viral_sig.xlsx')
        dfF1= pd.read_excel('Vi-Fi Signatures\\'+drug1+'_Fibrosis_sig.xlsx')
        dfF2= pd.read_excel('Vi-Fi Signatures\\'+drug2+'_Fibrosis_sig.xlsx')
        up_genes1= set(dict.fromkeys(dfV1.Gene[dfV1['change']=='up']))
        down_genes1= set(dict.fromkeys(dfV1.Gene[dfV1['change']=='down']))
        up_genes2= set(dict.fromkeys(dfV2.Gene[dfV2['change']=='up']))
        down_genes2= set(dict.fromkeys(dfV2.Gene[dfV2['change']=='down']))
        dfV1.set_index('Gene', inplace=True)
        dfV2.set_index('Gene', inplace=True)
        dfF1.set_index('Gene', inplace=True)
        dfF2.set_index('Gene', inplace=True)
        FS1= Score(dfF1, 'FScore')
        FS2= Score(dfF2, 'FScore')
        VS1= Score(dfV1, 'VScore')
        VS2= Score(dfV2, 'VScore')
        comb_up= list((up_genes1 | up_genes2) - (down_genes1 | down_genes2))
        comb_down= list((down_genes1 | down_genes2)-(up_genes1 | up_genes2))
        df_Vcomb=pd.concat([dfV1.loc[comb_up].dropna(), dfV2.loc[comb_up].dropna(), dfV1.loc[comb_down].dropna(), dfV2.loc[comb_down].dropna()])
        #df_comb=pd.concat([df1.loc[comb_up].dropna(), df2.loc[comb_up].dropna()])
        df_Vcomb.drop('Unnamed: 0', axis=1, inplace=True)
        df_Vcomb = df_Vcomb[~df_Vcomb.index.duplicated(keep='first')]
        df_Fcomb=pd.concat([dfF1.loc[comb_up].dropna(), dfF2.loc[comb_up].dropna(), dfF1.loc[comb_down].dropna(), dfF2.loc[comb_down].dropna()])
        df_Fcomb.drop('Unnamed: 0', axis=1, inplace=True)
        df_Fcomb = df_Fcomb[~df_Fcomb.index.duplicated(keep='first')]
        VS_comb= Score(df_Vcomb, 'VScore')
        FS_comb= Score(df_Fcomb, 'FScore')
        FS_list.append(FS_comb)
        FS1_list.append(FS1)
        FS2_list.append(FS2)
        VS_list.append(VS_comb)
        VS1_list.append(VS1)
        VS2_list.append(VS2)
        names_list.append(drug1+'&'+drug2)
    df_temp= pd.DataFrame({'Drugs combination': names_list, 
                               'FS_drug1': FS1_list, 'FS_drug2': FS2_list,
                               'FS_combination': FS_list, 'VS_drug1': VS1_list,
                               'VS_drug2': VS2_list, 'VS_combination': VS_list})
    df_score= pd.concat([df_score, df_temp])
df_score.to_excel(input_file+'_combinations.xlsx')

