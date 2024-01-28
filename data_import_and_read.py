# -*- coding: utf-8 -*-
"""
This module contains the functions used to import and translate data

Created on Fri Jan 28 18:24:41 2022

@author: Tacien Petithomme
@email: ptacien@gmail.com
License:  Creative Commons Attribution 4.0 International (CC BY 4.0)
Affiliation: INSERM 1307, Immunomodulation of the Tumor Microenvironment and Immunotherapy of Thoracic Cancers, Nantes, FRANCE




"""
import base64
import pandas as pd
import io
import dash_core_components as dcc
import dash_html_components as html

def read_file(contents, filename, remove_stop):
    """
    Reads the Dash-imported files and translates them into pandas dataframe. 
    The dataframe is then stored into a parquet file for convienience and fast accessibility in other process.
    
    Parameters
    ______________
    contents and filename are given by the Dash method for file import. 
    Notes that the function is able to read csv, xls, xlsx, fasta or txt files. In csv, xls, xlsx files, the file
    must contain only two columns. the fisrt must be the name of each sequence, and the second is the sequence in
    one string. 
    remove_stop is 0 or 1 from a Dash Checkbox. Informe whether or not to keep only sequences with no internal stop.
    
    Returns
    _____________
    Returns the Dataframe used in Dash to show the tables. 
    Stores the usefull dataframes in parquet files
    
    """
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
        # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
        elif 'xlsx' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.StringIO(decoded))
        elif 'fasta' in filename or 'txt' in filename:
            file=decoded.decode('utf-8')
            file=file.replace('\r', '')#sometimes, '\r' produces error
            l=file.split('>') #Split as a list of genes
            l = list(filter(None, l)) #this part is to remove empty object in the fasta file
            
            #then the list neads some clearing because of several issues that may be present the fasta files
            
            #sometimes, the sequence starts with a \n,  so i remove this artefact
            for i in range(len(l)):
                if l[i][0:2]=='\n': 
                    l[i]=l[i].replace('\n', '',1)
            #Now i'm supposed to have a string that contains the name and the sequence. 
                if '\n' in l[i]:
                    #I fisrt split the name from the sequence
                    l[i]=l[i].split('\n', maxsplit=1) 
                    #then I remove all \n inside the sequence
                    l[i][1]=l[i][1].replace('\n', '')
                    #if transform_stop==True:
                    #    l[i][1]=l[i][1].replace('*', 'Q')
            df=pd.DataFrame(l)
        
        #at this point, regardless the original file, i have a pandas dataframe
        df.columns=(['name', 'sequence'])
        #this part is to find sequences that carry ambigousn amino-acid caracters and remove them
        df=df[~df['sequence'].str.contains('X|J|B|Z')]
        
        #then we have to manage stop-containing seuqneces
        #to find and label stop-containing sequences:
        df['stop']=df.sequence.map(lambda x: '*' in x)
        
        if remove_stop==1: #alias if we want to remove them
            df=df.loc[df['stop']==False]
            df=df.drop(columns='stop')
        #but in bio-panning experiment, codon randomisation often use amber stop-codon as Gln (SupE/GlnV bacteria strain)
        else:
            df.sequence=df.sequence.map(lambda x: x.replace('*', 'Q'))
        df['file']=filename
        #know i store the parquet file.
        df.to_parquet('df.'+filename,engine='fastparquet', 
              compression='gzip') #write dataframe on local folder
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])
    return df


'''
the following function allows to adapt the Dash layout to the number of selected files
'''
def group_choice(contents, filenames):
    c=[]
    for i in filenames:
        c.append({'label':i, 'value':i})
    return dcc.Checklist(id='groups', options=c, value=filenames)


'''
This function take a dataframe and concatenates the rows so that each of them contains a UNIQUE sequence

def concatenate_unique(df):
'''



