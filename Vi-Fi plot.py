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
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

"""define working directory"""
os.chdir('C:\\Vi-Fi scoring')
df= pd.read_excel('Combined ViFi scoring.xlsx')
labels=list(dict.fromkeys(df.drug_class))
custom_color=['#b03315', '#e58725', '#FBC570','#ea96a3', '#ab9e47', '#97b9e3', '#7fa946', '#926ea9','#49d1b1', '#301D5B']
fig, ax = plt.subplots(figsize=(10, 10))
df1= df[(df['method']!='assay') & (df['method']!='mining')]
df2= df[df['method']=='assay']
df3= df[df['method']=='mining']
for dtype in list(dict.fromkeys(df1.drug_class_n)):
    df_plot=df1[df1['drug_class_n']==dtype]
    ax.scatter(df_plot['Viral Score'], df_plot['Fibrotic Score'], c=custom_color[dtype], s=150, marker='*', edgecolors='k')
for dtype in list(dict.fromkeys(df2.drug_class_n)):
    df_plot=df2[df2['drug_class_n']==dtype]
    ax.scatter(df_plot['Viral Score'], df_plot['Fibrotic Score'], c=custom_color[dtype], s=100, marker='o', edgecolors='k')
for dtype in list(dict.fromkeys(df3.drug_class_n)):
    df_plot=df3[df3['drug_class_n']==dtype]
    ax.scatter(df_plot['Viral Score'], df_plot['Fibrotic Score'], c=custom_color[dtype], s=100, marker='v', edgecolors='k')
#create legends
legend_elements1=[]
for dtype in list(dict.fromkeys(df.drug_class_n)):
    #legend_elements1.append(Circle((0,0), 1, color=custom_color[dtype], label=labels[dtype]))
    legend_elements1.append(Line2D([], [], marker='o', color=custom_color[dtype], linestyle='None', label=labels[dtype], markersize=10, markeredgecolor='k'))
legend1=ax.legend(handles=legend_elements1, loc='upper right', bbox_to_anchor=(1.12, 1.25), title="Drug type", prop={'size': 13}, edgecolor='#323232')
legend1.get_title().set_fontsize('15')
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
ax.add_artist(legend1)
legend_elements2 = [Line2D([], [], marker='o', color='#323232', linestyle='None', label='kinase assay', markersize=10),
                   Line2D([], [], marker='v', color='#323232', linestyle='None', label='text mining', markersize=10),
                   Line2D([], [], marker='*', color='#323232', linestyle='None', label='both methods', markersize=14)]
legend2=ax.legend(handles=legend_elements2, loc='upper right', bbox_to_anchor=(1.12, 0.83), title="Discovery method", prop={'size': 13}, edgecolor='#323232')
legend2.get_title().set_fontsize('15')
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
ax.add_artist(legend2)
#redo axis 
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
# remove the ticks from the top and right edges
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2)
# set axis labels
ax.set_xlabel('viral score', fontsize=18, weight='bold')
ax.xaxis.set_label_coords(1.0, 0.52)
ax.set_ylabel('fibrosis score', rotation=90, fontsize=18, weight='bold')
ax.yaxis.set_label_coords(0.46, 1.0)
#annotate drugs
annotation={'gefitinib':(-0.08, -0.005),
       'dasatinib':(-0.09, 0.01),
       'digoxin':(-0.07, -0.03),
       'piroxicam':(-0.09, -0.03),
       'hydrocortisone':(-0.12, -0.05),
       'estradiol':(-0.03, -0.05),
       'hg-6-64-01':(-0.09, -0.04),
       'calcitriol':(0.02, -0.03),
       'dexamethasone':(0.02, 0.02),
       'doxorubicin':(-0.1, -0.03),
       'proscillaridin':(-0.07, -0.05),
       'chloroquine':(-0.1, -0.04),
       'capsaicin':(0.03, -0.05)}
df_xy= df[['Viral Score', 'Fibrotic Score','Name']]
df_xy.set_index('Name', inplace=True)
for drug in list(dict.fromkeys(annotation)):
    xy = list(df_xy.loc[drug])
    ax.annotate(drug, xy, xytext=(xy[0]+annotation[drug][0], xy[1]+annotation[drug][1]), size=16, textcoords='data',
                arrowprops=dict(facecolor='#323232', width=1, headwidth=7, shrink=0))
fig.savefig('2D plot Fig5 small v2.pdf', bbox_inches = 'tight')