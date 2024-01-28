# -*- coding: utf-8 -*-
"""
This module contains the functions used to prossess the dataframe.
The first is here to split df.sequence to one amino-acid by column
The second/third calculate a consensus sequence
The last translates the sequence

Created on Fri Jan 28 18:24:41 2022

@author: Tacien Petithomme
@email: ptacien@gmail.com
License:  Creative Commons Attribution 4.0 International (CC BY 4.0)
Affiliation: INSERM 1307, Immunomodulation of the Tumor Microenvironment and Immunotherapy of Thoracic Cancers, Nantes, FRANCE


"""
import numpy as np
import pandas as pd
from constants import *




def split_df(df):
    '''
    Parameters
    ----------
    df : dataframe that must contain at least one columns with named 'sequence'

    Returns
    -------
    dfss : a pandas dataframe
        the first columns represent each amino-acid positions. All the columns that were
        not the sequence in df are added as last columns
    '''
    s=df['sequence'].to_numpy()
    if None in s:
        s=s[s!=None]
    dfss=pd.DataFrame(np.array(list(map(list, s))), dtype=(object))
    for i in df.columns: #petit ajout pour garder toutes les infos dans les colonnes autres que 'sequence'
        if i!='sequence':
            dfss[i]=df[i] # ATTENTION, parmis les i, il doit y avoir 'name' pour la suite des évènements
    return dfss


def get_repartition(df):
    '''
    this function is used in get_consensus(). 
    
    Parameters
    ----------
    df : dataframe with the amino acids sequences. One columns per position. Obtained from split_df3()

    Returns
    -------
    c.fillna(0) : a pandas dataframe with positions as column and all possible amino-acid as index. 
    the function counts the animo-acid occurence for each position and completes the dataframe
    '''
    c=pd.DataFrame(index=(aa))
    for i in df.columns: #je lis par colone
        if type(i) is int:
            c=pd.concat([c, df[i].value_counts()], axis=1)
    return (c.fillna(0))


def get_consensus(df):
    '''
    in order to convert categorical data (amino-acids) into numerical data, the strategy here is to 
    calculate a consensus sequence and express the dataset as a distance to consensus sequence. 
    
    The consensus is defined as the sequence with the shorter total distance to the entire
    dataset. 
    Distances are calculated based on grantham matrix. This matrix has the advantage on tacking into
    acount amino-acids physico-chemical features to calculate a distance between them. This matrix is more
    relevant in our matter as we are here interested in comparing final proteins obtained from random sequences. 
    Thus we don't need to calculate a phylogenic distance in term of number of mutation needed to pass from
    one amino-acid to another. 
    
    Parameters
    ----------
    df : dataframe with the amino acids sequences. One columns per position. Obtained from split_df()

    Returns
    -------
    df: the given dataframe
    consensus: the consensus sequence as a dict(). The keys are the position numbers and the values are
    the consensual amino-acid for a given position. 
    
    '''
    
    c=get_repartition(df).fillna(0)
    #c=cure_non_aa(c, df)
    consensus={}
    for pos in c:
        if len(np.where(c[pos].to_numpy()>0)[0])==1:
            consensus[pos]=aa[np.where(c[pos].to_numpy()>0)[0][0]]
        else:
            best=[]
            list_aa=[]
            for i in np.where(c[pos].to_numpy()>0)[0]:
                list_aa.append(aa[i])
                arr=(c[pos].to_numpy()*grantham_df[aa[i]].to_numpy())
                boolarr=arr>0
                best.append(arr[boolarr].mean())
            consensus[pos]=list_aa[np.where(best==min(best))[0][0]]
    return consensus


def matrix_df(df, consensus):
    '''
    This function translate the categorical dataframe (from split_df3) to a numerical dataframe
    by calculating for every positions in all the sequences the distance to the consensus sequence.
    The grantham matrix is formated as a dict() (grantham_d2) so the function simply applies the
    dict to the entire dataset.
    
    Parameters
    ----------
    df : dataframe with the amino acids sequences. One columns per position. Obtained from split_df()
    consensus : consensus sequence as a dict, obtained from consensus5()

    Returns
    -------
    A matrix that reprensent the dataset but where each amino-acid has been replaced by its grantham distance
    to the consensus sequence.
    '''
    matrix=pd.DataFrame()
    for pos in df.columns:
        if type(pos) is int:
            matrix[pos]=df.loc[:,pos].map(grantham_d2[consensus[pos]])
            if consensus[pos]=='-': #if '-' is the consensus sequence, I don't want any other AA to take the value 215
                #As is would loose a lot of information
                d=pd.DataFrame(df[pos])
                d=d.loc[d[pos]!='-']
                c=get_consensus(d) #I calculate the consensus AA without '-'
                nd=grantham_d2[c[1][pos]]
                nd['-']=0 #define a new distance dict
                d['m']=df.loc[:,pos].map(nd)
                matrix[pos]=matrix[pos]+d['m']
                matrix[pos]=matrix[pos].fillna(0)
    return matrix

def dist_matrix(df):
    '''this function takes a df from split function and returns a distance matrix'''
    d=pd.DataFrame()
    for i in range(df.shape[0]):
        d[i]=matrix_df(df, df.iloc[i,:].to_dict()).sum(axis=1)
    return(d)



'''
this function takes a split dataframe and reconvert it to inital imported dataframe.
'''
def split_to_total(df):
    r=pd.DataFrame()
    d=df
    for i in df.columns:
        if type(i) is not int:
            r[i]=df[i]
            d=d.drop(columns=i)
    d=d.to_csv(index=False)
    d=d.replace('\r', '')
    d=d.replace(',', '')
    l=d.split('\n')
    l=l[1:] #l[0] contains column names.
    dfx=pd.DataFrame(l)
    print(r.columns)
    for i in r.columns:
        dfx[i]=r[i]
    return(dfx)
    


