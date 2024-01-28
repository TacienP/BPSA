# -*- coding: utf-8 -*-
"""
This module contains the functions used to process data in order to plot them
Created on Fri Jan 28 18:24:41 2022

@author: Tacien Petithomme
@email: ptacien@gmail.com
License:  Creative Commons Attribution 4.0 International (CC BY 4.0)
Affiliation: INSERM 1307, Immunomodulation of the Tumor Microenvironment and Immunotherapy of Thoracic Cancers, Nantes, FRANCE




"""
import pandas as pd
import plotly
import plotly.graph_objects as go
from constants import *

'''
takes a dataframe with at least a 'file' and 'label' column and counts the number
of sequences in a file in every cluster. 
The output is a dataframe with every file as row and cluster number as column
'''
def cluster_rep(df):
    selection=df['file'].unique()
    label=list(df['label'].unique())
    label.sort()
    df_rep=pd.DataFrame(index=selection, columns=label)
    for l in label:
        d=[]
        for s in selection:
            part=len(df.loc[df['label']==l].loc[df['file']==s])
            d.append(part)
        df_rep[l]=d
    return(df_rep)

'''
This function takes a cluster_rep dataframe and produces a plotly sanky plot figure. 
'''
def sankey_plot(df):
    selection=[str(x) for x in df.index]
    cluster=[str(x) for x in df.columns]
    labels=selection+cluster
    sources=[]
    targets=[]
    values=[]
    for s in range(len(selection)):
        for c in range(len(cluster)):
            sources.append(s)
            targets.append(c+len(selection))
            values.append(df.iloc[s, c])
    fig = go.Figure(data=[go.Sankey(
        node = dict(
          pad = 15,
          thickness = 20,
          line = dict(width = 0.5),
          label = labels,
        ),
        link = dict(
          source =sources, # indices correspond to labels, eg A1, A2, A2, B1, ...
          target =targets,
          value =values
      ))],
        layout = go.Layout(height=800)
        )

    fig.update_layout(title_text="Cluster Sankey Diagram", font_size=10)
    return(fig)  


'''
This function counts the frequency of every amino-acid on every position. It returns a dataframe with
a compatible format for logomaker package.
'''
def lm_df2(df): #5 fois plus rapide. Est actuellement utilisé. 
    lmdict={}
    n=df.shape[0]
    for i in range(len(df.columns)):
        count=df[i].value_counts()
        l=[]
        for u in aa:
            if u in count:
                l.append(count[u]/n)
            else:
                l.append(0)
        lmdict[i]=l
    dflm=pd.DataFrame.from_dict(lmdict, orient='index', columns=aa)
    return dflm

def look_for_number(s, sep):
    l=s.split(sep)
    a=''
    for i in l[1]:
        if i.isdigit():
            a+=i 
    return (int(a))


'''
This function takes a dataframe

def sankey_inter_file(df):
 '''   

