# -*- coding: utf-8 -*-
"""
This module is the main one that contains the layout representation and executes all the functions.

Created on Sat Jan 29 16:45:43 2022

@author: Tacien Petithomme
@email: ptacien@gmail.com
License:  Creative Commons Attribution 4.0 International (CC BY 4.0)
Affiliation: INSERM 1307, Immunomodulation of the Tumor Microenvironment and Immunotherapy of Thoracic Cancers, Nantes, FRANCE


"""

'''
import section
'''
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import plotly.express as px
import os
import shutil
import pandas as pd
import tkinter
from pandas import ExcelFile, ExcelWriter
import numpy as np
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import math
import plotly
import plotly.graph_objects as go
import plotly.io as pio
import plotly.figure_factory as ff
from plotly.offline import plot
from plotly.subplots import make_subplots
import base64
import datetime
import io
import dash_table
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
from mpl_toolkits.mplot3d import Axes3D
from sklearn.neighbors import DistanceMetric
from sklearn.cluster import DBSCAN
import pyarrow
import logomaker as lm
import fastparquet
import webbrowser

from data_import_and_read import *
from matrix_calculation import *
from constants import *
from cluster_plot import *

'''
Working directory adjustment and creation of output folder.
'''
#choose your working directory
direction=os.getcwd()
if os.path.exists('BPSA_output'):
    os.chdir(direction+'/BPSA_output')
else:
    os.makedirs('BPSA_output')
    os.chdir(direction+'/BPSA_output')


'''
______________________________________________
Here a fisrt section with the Dash layout
______________________________________________
The callbachs are bellow
'''

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets,  suppress_callback_exceptions=True)

#this is how the app looks like
app.layout= html.Div([
    html.H1(children='Welcome to data analyzer'), #to be nice
    #the dashbord is constructed on several Tabs:
    dcc.Tabs([
            #this first tab allows user to import data, check that it worked properly and set whether to keep
            #stop-containing sequences or not(which replaces them by 'Q'). 
            #At the bottom of the tab, the user can run the matrix calculation and decide on which files the analysis 
            #must be done.
            dcc.Tab(label='Data import', children=[
            dcc.Upload(
                id='upload-data',
                children=html.Div([
                    'Drag and Drop or ',
                    html.A('Select a .fasta Files')
                ]),
                style={
                    'width': '100%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
                },
                # Allow multiple files to be uploaded
                multiple=True
            ),
            dcc.RadioItems(
                id='show_table',
                options=[{'label': 'Display table', 'value': 1}, 
                         {'label': 'Hide table', 'value': 0}],
                labelStyle={'display': 'inline-block'}, 
                value=1
            ),
            dcc.RadioItems(
                id='manage_stop',
                options=[{'label': 'Remove STOP-containing sequences', 'value': 1}, 
                         {'label': 'Keep STOP- containing sequences', 'value': 0}],
                labelStyle={'display': 'inline-block'},
                value=1
            ),
            html.Div(id='output-data-table'),
            ''' This was section was designed to allow further analysis to compare the evolution the clusters between rounds (i.e. files) and is not ready yet
            html.H6('Selection order'),
            html.P('If the files represent different selection tour, please indicate their number. Note that 2 files can have the same number for paralelle strategies'),
            html.P('Example : file_a.fasta / file_b.fasta / file_c.fasta / file_d.fasta  --> 1233'),
            html.P('In this example, file_a and file_b are the 2 fasta output of NGS of succissive selection tour while c and d are 2 different third tour strategies'),
            html.Div(id='inter'), 
            dcc.Input(id='order', type='text', placeholder='Example : 1233'),
            html.Hr(),
            ''',
            dcc.RadioItems(id='analysis',
                           options=[{'label': 'perfom analysis', 'value': 1}, 
                                    {'label': 'wait', 'value': 0}],
                labelStyle={'display': 'inline-block'},
                value=0),
            html.Div(id='group'),])
            ,
            
            #This second tab is to project the whole dataset on a single plot.
            #It is mostly to have a global view on data but no analysis are actually done here. 
            #As the dataset has been translated to a numerical matrix, is can be represented as a heatmap
            #or projected in two or tree-dimentional plot.
            #One may also decide on this tab to keep all amino-acid positions for further calculation
            #or remove some of them to speed up calculations if needed. An automatic removal of positions is also 
            #proposed and explained. 
            dcc.Tab(label='Whole data analysis', children=[
                    html.Details(html.P('''After importing the sequences, the files will be converted in a dataframe. Every column represents a position.
                                  A consensus will be calculated as the sequence with the minimal distance to the whole dataframe according to the Grantham matrix.
                                  The dataframe will the be translated in a matrix (refered as 'M' in the following details). 
                                  Each amino acid will take the value representing its distance to the consensus sequence (according to the Grantham matrix), at this position.''')
                                  ),
                    html.Div(id='matrix'),
                    html.P('Plot type: '),
                    dcc.Dropdown(id='plot_choice', options=[{'label': 'HEATMAP', 'value': 'heatmap'},
                     {'label': 'PCA 2D', 'value': 'PCA2D'},
                     {'label': 'PCA 3D', 'value': 'PCA3D'},
                     {'label':'show distance matrix', 'value':'show_matrix'}]),
                    html.P('Color map:'),
                    dcc.Dropdown(id='color'),
                    html.Div(id='plot',style = {'height' : 700}),
                    html.Hr(),
                    html.H6('Remove unvariable positions'), 
                    html.P('Please select you option'),
                    dcc.RadioItems(id='meth_removal', 
                                   options=[{'label':'Keep all positions', 'value':'all'},
                                            {'label':'Specify positions to keep', 'value':'keep'},
                                            {'label':'Specify positions to remove', 'value':'remove'},
                                            {'label':'Automatic removal', 'value':'auto'}],
                                   value='all'),
                    html.P('Type here the positions (starting at 0) you want to keep or remove in this format : 0,24,25,30'),
                    html.P('For automatic removal, please give the threshold (see details)'),
                    html.Div(id='detail_auto'),
                    dcc.Input(id='pos_removal', type='text', debounce=True),
                    html.Div(id='new_matrix'),
                    html.Hr(),
                    html.Div(id='new_plot',style = {'height' : 700}),
                    ])
            , 
            
            #this tab allows to cluster the data. 
            #On the top, the user can choose several clustering parameters. 
            #Then, the plots present the results of the clustering.
            dcc.Tab(label='All clusters', children=[
                        html.H6('Perform cluster analysis?'),
                        dcc.RadioItems(id='find clusters',
                           options=[{'label': 'no', 'value': 0}, 
                                    {'label': 'yes', 'value': 1},
                                    ], value=0,
                           labelStyle={'display': 'inline-block'}),
                        
                        html.P('Please set the DBSCAN parameters : '),
                        html.P('https://scikit-learn.org/stable/modules/clustering.html#dbscan'),
                        html.P('1- epsilon value // 2- minimum cluster population // 3- variance retained in dimention reduction (before analysis)'),
                        dcc.Input(id="e_value", type="number", placeholder="epsilon value", value=5),
                        dcc.Input(id="m_value", type="number", placeholder="min_sample", value=50),
                        dcc.Input(id="var", type="number", placeholder="variance retained (0 - 1)", value=0.8),
                        dcc.Checklist(id='consider_number',
                                      options=[{'label': 'Take number of sequences into acount', 'value': 'yes'}]),
                        dcc.Input(id="separator", type="text", placeholder="separator in sequences name"),
                        html.Details(html.P('''Sometimes, the files contain the list of non-redundant sequences
                                            and the name of each sequence contains the number of times this sequence was found.
                                            This is a function to re-edit a matrix containing the full dataset.
                                            Separator is the specific string characters found just before the number information 
                                            in the name of each sequence. Ex : 
                                            sequence No1. ejch ##440##
                                            sequence No2. ejch ##30##
                                            A possible separator to indicate where is the number in the name is "##"''')),
                        html.Div(id='clusters_calculation'),
                        html.Div(id='num_of_clusters'),
                        html.P('Display results?'),
                        dcc.RadioItems(id='display clusters',
                           options=[{'label': 'No', 'value': 0}, 
                                    {'label': 'Yes', 'value': 1},
                                    ], value=0,
                           labelStyle={'display': 'inline-block'}),
                        dcc.RadioItems(id='how_to_plot_cluster',
                           options=[{'label': 'PCA', 'value': 0}, 
                                    {'label': 'T-SNE', 'value': 1}, 
                                    {'label': 'UMAP', 'value': 2}
                                    ], value=0,
                           labelStyle={'display': 'inline-block'}),
                        dcc.RadioItems(id='plot_cluster_dimention',
                           options=[{'label': '2D', 'value': 2}, 
                                    {'label': '3D', 'value': 3}, 
                                    ], value=2,
                           labelStyle={'display': 'inline-block'}),
                        dcc.Checklist(id='plot_only_cluster',
                                      options=[{'label': 'Plot only clustered sequences', 'value': 'yes'}],
                                      value=['yes']),
                        html.Div(id='plot_clusters',style = {'height' : 1000}),
                        dcc.Graph(id='scatter_contour'),
                        html.Div(id='clusters_heatmap'),
                    ])
            ,
            
            #This tab is specialized to analyse the variation between different files. 
            #Note that file must correspond to different selection rounds.
            #Selection order must be mentionned on the fisrt tab. 
            dcc.Tab(label='Inter-files analysis', children=[
                        html.P('Read repartition between files. Note that equal repartition is needed to analyze the evolution of the populations between rounds'),
                        html.Div(id='reads_count'),
                        html.Div(id='bar_clusters'),
                    
                    ])
            , 
            
            #In this tab, the user can specify the cluster to be displayed.
            dcc.Tab(label='Intra cluster', children=[
                        html.P('You can select a region on this plot to retrieve there sequences'),
                        dcc.Graph(id='scatter_contour2'),
                        html.Div(id='selected_points_table'),
                        html.P('To analyse a particular cluster, please enter the cluster number here:'),
                        dcc.Input(id='selecte_cluster', type='number'), #select a cluster
                        html.Div(id='cluster_table'), #exportable table of the sequences in the selected cluster
                        html.Div(id='cluster_pie'),#file repartition in this cluster
                        html.Div(id='cluster_plot'), #re-analayze and plot the sequences in this cluster
                         ])
            ,
            ], vertical=False)
    ])
                                            

'''
______________________________________________
Second section with the Dash callbacks
______________________________________________
'''

#callback to manage uploaded dataframe. It also displays the dashTable
@app.callback(Output('output-data-table', 'children'),
              Input('upload-data', 'contents'),
              Input('show_table', 'value'),
              Input('manage_stop', 'value'),
              State('upload-data', 'filename'))
def update_output(list_of_contents, display_option, stop, list_of_names ):
    if list_of_contents is None:
        return html.H6('Please load data')
    else:
        #read files and transforme them in pandas dataframes
        #in the loop, read_file stores each file's dataframe in a parquet
        all_df = [read_file(c, n, stop) for c, n in zip(list_of_contents, list_of_names)]
        #the parquet are then retrieved to creat a concatenated dataframe (total_df)
        total_df=pd.DataFrame()
        for i in list_of_names:
            dft=pd.read_parquet('df.'+i)
            total_df=pd.concat([dft, total_df], ignore_index=True)
        #total_df is stored in parquet format
        total_df.to_parquet('total_df',engine='fastparquet', 
                  compression='gzip')
        if display_option==1:
            res=[
                html.Div([
                html.H5(filename),
                    dash_table.DataTable(
                        id='table',
                        data=(data.to_dict('records')),
                        columns=[{'name': i, 'id': i} for i in data.columns],
                        page_size=10,
                        export_format="csv",
                        ),
                html.Hr(),  # horizontal line
                        ])
                for filename, data in zip(list_of_names, all_df)]
            
            return res
        else:
            return html.H6('Data uploaded')




#This callback allows to display the name of the different files. 
@app.callback(Output('inter', 'children'), 
              Input('upload-data', 'filename'))
def get_order(list_of_names):
    if list_of_names is not None:
        return html.P(" | ".join((str(name) for name in list_of_names)))


#this callback generates a checkbox that allows the user to select the files that have to be analysed.
@app.callback(Output('group', 'children'),
              Input('upload-data', 'contents'),
              Input('analysis', 'value'),
              State('upload-data', 'filename'),
              suppress_callback_exceptions=True)
def check_menu(list_of_contents, perform, list_of_names):
    if perform==1:
        return group_choice(list_of_contents, list_of_names)


#this callbak is used to compute matrix once for all and store it on the computer
@app.callback(Output('matrix', 'children'),
              Input('groups', 'value'),
              prevent_initial_call=True,
              suppress_callback_exceptions=True)
def generate_matrix(group):
    '''
    In this funstion, the categorical data are translated into numerical data. 
    Note that the output matrix is not the one that is used eslewhere in the programme
    since the following callback (new_matrix()) will potentially revome positions with low variability. 
    '''
    if group is not None:
        dfa=pd.read_parquet('total_df')
        #This part keeps the files selected in check_menu.
        analyse_df=pd.DataFrame()
        u=[i for i in group]
        analyse_df=dfa[[x in u for x in dfa.file]]
        #Store the data that actually have to be analysed
        analyse_df.to_parquet('pre-matrix.gzip', engine='fastparquet', compression='gzip')
        #_______________________________________________
        #in case stop codon were considered as Gln, a column named 'stop' 
        #is in total_df. Lets first remove it in matrix.
        if 'stop' in dfa:
            analyse_df=analyse_df.drop(columns=['stop'])
        #_______________________________________________
        #This part translates the data into a numeric matrix
        split=split_df(analyse_df)
        consensus=get_consensus(split)
        matrix=matrix_df(split, consensus)
        matrix['file']=analyse_df['file']
        #save splited and matrix dataframe on local folder for further use
        split.columns=[str(j) for j in list(split.columns)]
        split.to_parquet('split_df.gzip', engine='fastparquet', compression='gzip')
        matrix.columns=[str(j) for j in list(matrix.columns)]
        matrix.to_parquet('matrix0.gzip', engine='fastparquet', compression='gzip')
        
        #A little part to let the user know what is happening.
        if None in analyse_df.sequence:
            a='Be carreful, "None" sequece found in the file and it may induce frame decay. Please check the files.'
        else:
            a='Calculations done. You may plot them below.'
        #The only return of this function is a. Important calculation are stored in parquet at this stage.
        return html.H6(a)


#callback that presents the representation options
@app.callback(Output('color', 'options'),
              Input('manage_stop', 'value'),
              Input('plot_choice', 'value'),
              Input('color', 'value')
              )
def plot_color(stop, graph, cur):
    options=[]
    if graph=='PCA2D' or graph=='PCA3D':
        if stop==0:
            options.append({'label':'Stop-containing sequences', 'value':'stop'})
        options.append({'label':'By files', 'value':'files'})
    for i in px.colors.named_colorscales():
        options.append({'label':i, 'value':i})
    return options

#This function offers the user to visualize the variability of each position and explains how the automatic seleciton
#of variable position works. It helps the user to remove unvariable positions.
@app.callback(Output('detail_auto', 'children'),
              Input('matrix', 'children'),
              Input('upload-data', 'filename'))
def analyze_variability(m, imported):
    if m is not None:
        matrix=pd.read_parquet('matrix0.gzip')
        matrix=matrix.iloc[:,:-1] #remove 'file' column
        l=[]
        for i in matrix:
            s=matrix[i].sum()
            l.append(s)
        tot=sum(l)
        pie_var=px.pie(names=matrix.columns, values=l, title='Pie chart of sequences variability, by position')
        bar_var=px.bar(y=l, title='Bar chart of sequences variability, by position')
        ''' 
        nice but very long to show
        coll=l
        coll=np.array(coll).reshape(-1, 1)
        scaler = MinMaxScaler()
        coll=scaler.fit_transform(coll)
        colors=['rgb('+str(10+int(i))+','+str(220-int(i))+','+str(220-int(i))+')' for i in 200*coll]
        violin_var=go.Figure()
        for data_line, col in zip(matrix.columns,colors):
            violin_var.add_trace(go.Violin(y=matrix[data_line], line_color=col, name=data_line))
        violin_var.update_layout(yaxis_title='grantham matrix value')
        '''
        return(html.Div([
            html.Details([
                html.Hr(),
                html.H6('''the variability is measured here by calculating the sum of every column of the matrix M.
                   Total variability is the sum of every cells over the matrix.'''),
                dcc.Graph(figure=pie_var),
                html.Hr(),
                dcc.Graph(figure=bar_var),
                #dcc.Graph(figure=violin_var),
                html.P('''The threshold  used below is the percentage of total variability to 
                   reach for a position to be considered as a variable enough position. Typing '1' deletes all columns of M
                   with sum(column)<0.01*(total_sum(M))'''),
                         ]),
            ])
                   )


'''
This functions modifies the matrix (stored as matrix0.gzip) and only keeps the positions decided by the user. 
In case all the positions are kept for analysis, the new_matrix is simply the same as the one 
already computed. 
At the end of this function, a 3-dimention PCA reduction is computed once for all and stored 
so that ploting further function will only have to read the stored parquet.
'''
@app.callback(Output('new_matrix', 'children'),
              Input('pos_removal', 'value'), #text input.
              Input('meth_removal', 'value'), #tells whether to keep or remove.
              Input('matrix', 'children')) #allows to execute only after the first matrix has alread been computed.
def new_matrix(where,how, done):
    if done is not None:
        
        
        
        #By default, all the positions are kept, so there is no modifications to do on matrix0.gzip
        if how=='all':
            shutil.copy2('matrix0.gzip', 'matrix.gzip') #if no modification -> matrix=matrix0
            X=pd.read_parquet('matrix0.gzip').iloc[:, :-1] #the last column is the file name so we remove it.
        #______________________________________________
        #If the user wants to remove some amino-acids positions, there is something in "pos_removal"
        elif where is not None: #in this case, we will have to work on matrix0.gzip so lets start!
            m=pd.DataFrame() # adataframe that will contain only selected data
            matrix=pd.read_parquet('matrix0.gzip')
            #_ _ _ _ _ _ _ _ _ _ _ _ _ _
            if how=='keep':
                pos=where.split(',') #list of interesting columns
                for i in pos:
                    m[i]=matrix[i]
                m['file']=matrix['file']
            elif how=='remove':
                pos=where.split(',') #list of undesired data
                col=matrix.columns
                for i in col:
                    if i not in pos:
                        m[i]=matrix[i]
            elif how=='auto':
                l={}
                for i in matrix.iloc[:,:-1]:
                    s=matrix[i].sum()
                    l[i]=s
                tot=sum(list(l.values()))
                for i in l.keys():
                    if l[i]>(((float(where))/100)*tot): #if variability reaches the threshold, it is copied in m.
                        m[i]=matrix[i]
                m['file']=matrix['file']
            #At this point, we can save the new matrix.
            m.to_parquet('matrix.gzip', engine='fastparquet', compression='gzip')
            X=m.iloc[:, :-1] #removes file column
        #______________________________________________
        # Now we will compute PCA reduction on cleaned data once and for all
        pca_plot = PCA(n_components=3)
        PC_plot3D = pd.DataFrame(pca_plot.fit_transform(X))
        PC_plot3D.columns=['PC1', 'PC2', 'PC3']
        PC_plot3D.to_parquet('PC_plot_3D.qzip', engine='fastparquet', compression='gzip')
        return()


'''
We are in the 'Whole data analysis' tab
This callback allow to plot the dataset after matrix0.gzip it has been generated.
A first section imports PCA-reduced data then the plot are generated. 
'''
@app.callback(Output('plot', 'children'),
              Input('plot_choice', 'value'),
              Input('matrix', 'children'),
              Input('color', 'value'),
              )
def plot_data(plot_type, analysis, colorv):
    if analysis is not None:
        if plot_type in ['heatmap', 'PCA2D', 'PCA3D', 'show_matrix']:
            matrix=pd.read_parquet('matrix0.gzip')
            X=matrix.iloc[:, :-1]
            #___________HEATMAP PLOTTING____________
            if plot_type=='heatmap':
                if colorv=='files':
                    col='viridis'
                else:
                    col=colorv
                fig = go.Figure(data=go.Heatmap(z=X, colorscale=col), 
                                layout = go.Layout(height=700),
                                )
                return dcc.Graph(figure=fig)
            
            #___________2D PCA PLOTTING____________
            elif plot_type=='PCA2D':
                #First, read data
                df_to_plot=pd.read_parquet('PC_plot_3D.qzip')
                #Then, manage color assignment
                if colorv not in ['files', 'stop']:
                    m=list(range(matrix.shape[0]))
                    colscale=colorv
                elif colorv=='files':
                    col=[i for i in range(len(matrix['file'].unique()))]
                    c={i:j for i, j in zip(matrix['file'].unique(), col)}
                    m=matrix['file'].map(c) #file names are replaced by numbers.
                    colscale='Rainbow'
                elif colorv=='stop':
                    stop_df=pd.read_parquet('pre-matrix.gzip')['stop']
                    c={False:0, True:1}
                    m=stop.map(c)
                #Finally, plot!    
                fig=go.Figure(data=go.Scattergl(x=df_to_plot['PC1'], y=df_to_plot['PC2'], 
                                                mode='markers', 
                                                marker=dict(color=m,
                                                            colorscale=colscale, size=2)
                                                )
                              ,layout = go.Layout(height=700)
                              )
                return dcc.Graph(figure=fig)
                
            #___________3D PCA PLOTTING____________
            elif plot_type=='PCA3D':
                #First, read data
                df_to_plot=pd.read_parquet('PC_plot_3D.qzip')
                #Then, manage color assignment
                if colorv not in ['files', 'stop']:
                    m=list(range(matrix.shape[0]))
                    colscale=colorv
                elif colorv=='files':
                    col=[i for i in range(len(matrix['file'].unique()))]
                    c={i:j for i, j in zip(matrix['file'].unique(), col)}
                    m=matrix['file'].map(c) #file number are replaced by numbers.
                    colscale='Rainbow'
                elif colorv=='stop':
                    stop_df=pd.read_parquet('pre-matrix.gzip')['stop']
                    c={False:0, True:1}
                    m=stop.map(c)
                #Finally, plot!  
                fig=go.Figure(data=[go.Scatter3d(x=df_to_plot['PC1'], y=df_to_plot['PC2'],
                                                 z=df_to_plot['PC3'], mode='markers', 
                                                 marker=dict(size=2,color=m, 
                                                             colorscale=colorv,opacity=0.8)
                                                 )
                                    ]
                              ,layout = go.Layout(height=700)
                              )
                return dcc.Graph(figure=fig)
            
            else:
                matrix=pd.read_parquet('matrix.gzip')
                return dash_table.DataTable(
                    id='table2',
                    data=matrix.to_dict('records'),
                    columns=[{'name': i, 'id': i} for i in [str(j) for j in list(matrix.columns)]],
                    page_size=10
                    )
            
'''
We are in the 'Whole data analysis' tab
This callback allow to plot the dataset after the new matrix.gzip it has been generated.
A first section imports PCA-reduced data then the plot are generated. 
'''
@app.callback(Output('new_plot', 'children'),
              Input('plot_choice', 'value'),
              Input('matrix', 'children'),
              Input('color', 'value'),
              Input('new_matrix', 'children')
              )
def plot_data2(plot_type, analysis, colorv, new):
    if analysis is not None:
        matrix=pd.read_parquet('matrix.gzip')
        X=matrix.iloc[:, :-1]
        #___________HEATMAP PLOTTING____________
        if plot_type=='heatmap':
            if colorv=='files':
                col='viridis'
            else:
                col=colorv
            fig = go.Figure(data=go.Heatmap(z=X, colorscale=col), 
                            layout = go.Layout(height=700),
                            )
            return dcc.Graph(figure=fig)
        
        #___________2D PCA PLOTTING____________
        elif plot_type=='PCA2D':
            #First, read data
            df_to_plot=pd.read_parquet('PC_plot_3D.qzip')
            #Then, manage color assignment
            if colorv not in ['files', 'stop']:
                m=list(range(matrix.shape[0]))
                colscale=colorv
            elif colorv=='files':
                col=[i for i in range(len(matrix['file'].unique()))]
                c={i:j for i, j in zip(matrix['file'].unique(), col)}
                m=matrix['file'].map(c) #file names are replaced by numbers.
                colscale='Rainbow'
            elif colorv=='stop':
                stop_df=pd.read_parquet('pre-matrix.gzip')['stop']
                c={False:0, True:1}
                m=stop.map(c)
            #Finally, plot!    
            fig=go.Figure(data=go.Scattergl(x=df_to_plot['PC1'], y=df_to_plot['PC2'], 
                                            mode='markers', 
                                            marker=dict(color=m,
                                                        colorscale=colscale, size=2)
                                            )
                          ,layout = go.Layout(height=700)
                          )
            return dcc.Graph(figure=fig)
            
        #___________3D PCA PLOTTING____________
        elif plot_type=='PCA3D':
            #First, read data
            df_to_plot=pd.read_parquet('PC_plot_3D.qzip')
            #Then, manage color assignment
            if colorv not in ['files', 'stop']:
                m=list(range(matrix.shape[0]))
                colscale=colorv
            elif colorv=='files':
                col=[i for i in range(len(matrix['file'].unique()))]
                c={i:j for i, j in zip(matrix['file'].unique(), col)}
                m=matrix['file'].map(c) #file number are replaced by numbers.
                colscale='Rainbow'
            elif colorv=='stop':
                stop_df=pd.read_parquet('pre-matrix.gzip')['stop']
                c={False:0, True:1}
                m=stop.map(c)
            #Finally, plot!  
            fig=go.Figure(data=[go.Scatter3d(x=df_to_plot['PC1'], y=df_to_plot['PC2'],
                                             z=df_to_plot['PC3'], mode='markers', 
                                             marker=dict(size=2,color=m, 
                                                         colorscale=colorv,opacity=0.8)
                                             )
                                ]
                          ,layout = go.Layout(height=700)
                          )
            return dcc.Graph(figure=fig)
        
        else:
            matrix=pd.read_parquet('matrix.gzip')
            return dash_table.DataTable(
                id='table2',
                data=matrix.to_dict('records'),
                columns=[{'name': i, 'id': i} for i in [str(j) for j in list(matrix.columns)]],
                page_size=10
                )

'''
We are now in the cluster analysis tab
this callabck allows to run sequence clusterisation via DBSCAN method. 
'''
@app.callback(Output('clusters_calculation', 'children'),
              Input('plot_cluster_dimention', 'value'), 
              Input('e_value', 'value'), 
              Input('m_value', 'value'),
              Input('var', 'value'),
              Input('analysis', 'value'), 
              Input('matrix', 'children'),
              Input('find clusters', 'value'),
              Input('consider_number', 'value'),
              Input('separator', 'value'),
              )
def find_clusters(dim, e, m, v, analysis, done, find, num, sep):
    if analysis==1 and find==1:
        # Import data
        matrix=pd.read_parquet('matrix.gzip')
        x=matrix.iloc[:, :-1]
        # Run PCA reduction with desired retained variability. This allows to speed up the DBSCAN
        pca_cluster = PCA(v)
        PC_cluster = pca_cluster.fit_transform(x)
        ''' The clustering can be done on unique sequences or NGS reads.  
        Imported data may be the list of read or a concatenated version in which the number of time
        each unique sequence was observed in the dataset is included in the name of the sequence. 
        Annotation may vary but for example:
        >sequence no.1__##230
        ADKDEPRMEKDADDDFHETPSNCHKMATYUICNF
        The number after '##' indicates that this sequence was observed 230 times. 
        The web app allows the user to indicate the string sequence preceding the number.
        '''
        # If number of time each sequence are present are to be considered:
        if num==['yes']:
            if sep is not None:
                dft=pd.read_parquet('pre-matrix.gzip')
                PC_cluster=np.append(PC_cluster, np.array(dft.name).reshape(len(dft.name), 1), axis=1)
                PC_cluster=np.append(PC_cluster, np.array(dft.file).reshape(len(dft.file), 1), axis=1)
                try:
                    dft['number']=dft.name.map(lambda x: look_for_number(x, sep)) 
                except:
                    dft.name=dft.name.map(lambda x: str(x)+sep+'200' if sep not in x else x)
                    dft['number']=dft.name.map(lambda x: look_for_number(x, sep))
                    print("there was an error! Be carreful, every sequence name must contain the separator. A default '200' value was attributed to missing data")
                PC_cluster=PC_cluster.repeat(dft.number, axis=0)
                db=DBSCAN(eps=e, min_samples=m).fit(PC_cluster[:, :-2])
                PC_cluster=np.append(PC_cluster, np.array(db.labels_).reshape(len(db.labels_), 1), axis=1)
                d=pd.DataFrame(PC_cluster[:, -3:]).drop_duplicates()
                #a='Calculation done. I took the number into acount'
                for i in range(len(d[0])):
                    if d.iloc[i, 0]!=dft.name[i]:
                        a='Warning, label number and sequence name are not correctly aligned'
                        '''you should therefor manually use the d dataframe where name and label are aligned'''
                df_cluster=pd.DataFrame(d[2])
                df_cluster.columns=['label']
                df_cluster.index=range(df_cluster.shape[0])
                a='Reads were clustered in {} clusters'.format((len(df_cluster['label'].unique())-1))
            else:
                a='Please indicate a separator (see details)'

        # If we only want to cluster unique sequences
        else: #thereby num is not checked
            db=DBSCAN(eps=e, min_samples=m).fit(PC_cluster)
            df_cluster=pd.DataFrame(db.labels_, columns=['label'])
            a='Unique sequences were clustered in {} clusters'.format((len(df_cluster['label'].unique())-1))
        df_cluster.to_parquet('df_cluster.gzip',engine='fastparquet', 
                  compression='gzip')
        
        return html.H6(a)
         

'''
@app.callback(Output('num_of_clusters', 'children'),
              Input('analysis', 'value'), 
              Input('clusters_calculation', 'children'),
              Input('consider_number', 'value'),
              )
def count_cluster(analysis, done, num):
    if analysis==1 and done is not None:
        df_labels=pd.read_parquet('df_cluster.gzip')
        label=df_labels['label'].unique()
        if num is not None:
            a='Reads were clustered in {} clusters'.format((len(label)-1))
        elif num=='yes':
            a='Unique sequences were clustered in {} clusters'.format((len(label)-1))
        return html.H6(a)
    '''

'''
This callback plots the clustered sequences as a scatter plot. Note that only sequences that were 
regrouped as clusters will be represented. 
'''
@app.callback(Output('plot_clusters', 'children'),
              Input('how_to_plot_cluster', 'value'), #0=PCA, 1=T-SNE, 2=UMAP
              Input('plot_cluster_dimention', 'value'),
              Input('var', 'value'),
              Input('analysis', 'value'),
              Input('clusters_calculation', 'children'),
              Input('display clusters', 'value'),
              Input('plot_only_cluster', 'value'),
              )
def plot_cluster(plot_type, dim, v, analysis, done, display, plot_only_cluster):
    if analysis==1 and done is not None and display==1:
        matrix=pd.read_parquet('matrix.gzip')
        x=matrix.iloc[:, :-1]
        df_labels=pd.read_parquet('df_cluster.gzip')
        if plot_type==0:
            df_to_plot=pd.read_parquet('PC_plot_3D.qzip')
        elif plot_type==1: #for T-SNE and UMAP, to accelerate the calulation, we will performe the projection after a PCA reduction
            pca_cluster = PCA(v)
            PC_cluster = pca_cluster.fit_transform(x)
            plot_reduction=TSNE(n_components=dim)
            PC_plot=plot_reduction.fit_transform(PC_cluster)
            df_to_plot = pd.DataFrame(data = PC_plot
                     , columns = ['PC{}'.format(i) for i in range(dim)])
        elif plot_type==2:
            pca_cluster = PCA(v)
            PC_cluster = pca_cluster.fit_transform(x)
            plot_reduction=umap.UMAP(n_components=dim)
            PC_plot=plot_reduction.fit_transform(PC_cluster)
            df_to_plot = pd.DataFrame(data = PC_plot
                     , columns = ['PC{}'.format(i+1) for i in range(dim)])
        
        df_to_plot['label']=df_labels['label']
        if dim==2:
            fig=go.Figure(data=go.Scattergl(x=df_to_plot.loc[df_to_plot['label']!=-1]['PC1'],
                                            y=df_to_plot.loc[df_to_plot['label']!=-1]['PC2'], 
                                            mode='markers', 
                                            marker=dict(color=df_to_plot.loc[df_to_plot['label']!=-1]['label'],
                                                        colorscale='Viridis',size=3),
                                            text=df_to_plot.loc[df_to_plot['label']!=-1]['label'])
                          ,layout = go.Layout(height=1000))
            if plot_only_cluster!=['yes']:
                fig.add_trace(go.Scattergl(x=df_to_plot.loc[df_to_plot['label']==-1]['PC1'],
                                           y=df_to_plot.loc[df_to_plot['label']==-1]['PC2'], 
                                           mode='markers', 
                                           marker=dict(color='Grey',size=3,opacity=0.1),
                                           text='Unclustered data'
                                           )
                              )
            
        elif dim==3:
            # Here we plot clustered data
            fig=go.Figure(data=go.Scatter3d(x=df_to_plot.loc[df_to_plot['label']!=-1]['PC1'],
                                            y=df_to_plot.loc[df_to_plot['label']!=-1]['PC2'], 
                                            z=df_to_plot.loc[df_to_plot['label']!=-1]['PC3'],
                                            mode='markers', 
                                            marker=dict(color=df_to_plot.loc[df_to_plot['label']!=-1]['label'],
                                                        colorscale='Viridis',size=3,opacity=0.8),
                                            text=df_to_plot.loc[df_to_plot['label']!=-1]['label'])
                          ,layout = go.Layout(height=1000))
            # And then we may add unclustered data
            if plot_only_cluster!=['yes']:
                fig.add_trace(go.Scatter3d(x=df_to_plot.loc[df_to_plot['label']==-1]['PC1'],
                                           y=df_to_plot.loc[df_to_plot['label']==-1]['PC2'], 
                                           z=df_to_plot.loc[df_to_plot['label']==-1]['PC3'],
                                           mode='markers', 
                                           marker=dict(color='Grey',size=1,opacity=0.1),
                                           text='Unclustered data'
                                           )
                              )
        fig.update_layout(template='plotly_white')
        return dcc.Graph(figure=fig)

'''
This function generates several plots that intend to give a good representation and comprehension 
of the data and the different families(=clusters) of sequences
'''
@app.callback(Output('clusters_heatmap', 'children'),
              Input('analysis', 'value'),
              Input('clusters_calculation', 'children'),
              Input('display clusters', 'value'),
              )
def plot_cluster(analysis, done, display):
    if analysis==1 and done is not None and display==1:
        #Here a first part to generate a heatmap reprensetation of the the consensus
        # sequence of every cluster. 
        df_labels=pd.read_parquet('df_cluster.gzip') #to get cluster affiliation
        dfl=pd.read_parquet('split_df.gzip') # We go back to the splited dataframe
        col=[]
        for i in dfl.columns: #this is because column name are str instead of int after parquet saving...
            try :
                col.append(int(i))
            except:
                col.append(i)
        dfl.columns=col
        df_all=dfl.iloc[:,:-2] #removing 'file' and 'name', allows to save row data for later use. 
        dfl['label']=df_labels['label']
        df_count=dfl.loc[:,['file', 'label']] #...seems useless
        dfl=dfl.loc[dfl['label']!=-1] # We select only clustered sequences
        lbl=list(dfl.label.unique())
        lbl.sort()
        rows={}
        for i in lbl:
            rows[i]=list(get_consensus(dfl.loc[dfl['label']==i]).values()) #gets the consensus sequence of every cluster
        df_cc=pd.DataFrame.from_dict(rows, orient='index') #this dataframe contains the consensus sequence of every cluster
        #Now we generate a matrix from df_cc...
        consensus=get_consensus(df_cc)
        matrix=matrix_df(df_cc, consensus)
        #...and plot a heatmap
        fig_heatmap = go.Figure(data=go.Heatmap(z=matrix, colorscale='Viridis',
                                        customdata = df_cc,
                                        hovertemplate='<br>amino acid: %{customdata} '), 
                                        layout = go.Layout(height=800),
                                        )
        #___________________________________
        #Generate distance matrix between every clusters
        d_matrix=dist_matrix(df_cc)
        fig_dist = go.Figure(data=go.Heatmap(z=d_matrix, x=df_cc.index, y=df_cc.index))
        
        #___________________________________
        #generate pie chart
        '''
        l=list(df_labels['label'].unique())
        l.sort()
        lc=[i for i in l]
        '''
        lc=lbl.copy()
        lc.remove(-1)
        val=[len(df_labels.loc[df_labels['label']==i]), len(df_labels.loc[df_labels['label']!=i])]
        valc=[len(df_labels.loc[df_labels['label']==i]) for i in lc]
        p=[0.2, 0]
        pie_all=go.Figure(data=[go.Pie(labels=['Clusted', 'Background'], values=val, pull=p)],
                      layout = go.Layout(height=1000))
        pie=go.Figure(data=[go.Pie(labels=lc, values=valc)],
                      layout = go.Layout(height=1000))
        
        #____________________________________
        #generates sankey plot
        df_cluster_rep=cluster_rep(dfl)
        fig_sankey=sankey_plot(df_cluster_rep)
        #____________________________________
        #generates dendrogram
        dend=ff.create_dendrogram(d_matrix)
        dend.update_layout(height=600, width=1800)
        #____________________________________
        # generate fancy logo diagrams 
        plt.ioff()
        #compute and save whole dataset logo
        lm_all=lm_df2(df_all)
        all_im=lm.Logo(lm_all,color_scheme='dmslogo_funcgroup')
        all_im.fig.savefig('all_img.png',dpi=1200)
        #all_im_filename =direction+'\\all_img.png'
        all_im_encoded_image = base64.b64encode(open('all_img.png', 'rb').read())
        plt.close()
        
        #compute and save clustered sequence logo
        df_clustered=dfl.iloc[:,:-3] # remove 'file', 'name' and 'label' columns
        lm_clustered=lm_df2(df_clustered)
        cl_im=lm.Logo(lm_clustered,color_scheme='dmslogo_funcgroup')
        cl_im.fig.savefig('cl_img.png',dpi=1200)
        #cl_im_filename =direction+'\\cl_img.png'
        cl_im_encoded_image = base64.b64encode(open('cl_img.png', 'rb').read())
        plt.close()
        
        #compute and save cluster logo
        lm_cc=lm_df2(df_cc)
        cluster_im=lm.Logo(lm_cc,color_scheme='dmslogo_funcgroup')
        cluster_im.fig.savefig('cluster_img.png',dpi=1200)
        #cluster_im_filename =direction+'\\cluster_img.png'
        cluster_im_encoded_image = base64.b64encode(open('cluster_img.png', 'rb').read())
        plt.close()
        
        #___________________________________
        # Now, return the figures.
        return html.Div([
            html.Div([
                dcc.Graph(figure=pie_all, style={'display': 'inline-block'}),
                dcc.Graph(figure=pie, style={'display': 'inline-block'})
                ]),
            dcc.Graph(figure=fig_sankey),
            dcc.Graph(figure=fig_heatmap),
            dcc.Graph(figure=dend),
            html.H6('Whole dataset Logo'),
            html.Img(src='data:image/png;base64,{}'.format(all_im_encoded_image.decode()),  
                     style={'width': '100%'}),
            html.H6('Clustered sequences Logo'),
            html.Img(src='data:image/png;base64,{}'.format(cl_im_encoded_image.decode()),  
                     style={'width': '100%'}),
            html.H6("Cluster's consensus Logo"),
            html.Img(src='data:image/png;base64,{}'.format(cluster_im_encoded_image.decode()),  
                     style={'width': '100%'}),
            dash_table.DataTable(
                    id='table3',
                    data=df_cluster_rep.to_dict('records'),
                    columns=[{'name': i, 'id': i} for i in [str(j) for j in list(df_cluster_rep.columns)]],
                    page_size=10
                    )
            
            ])



@app.callback([Output('scatter_contour', 'figure'),
              Output('scatter_contour2', 'figure')],
              [Input('analysis', 'value'),
              Input('clusters_calculation', 'children'),
              Input('display clusters', 'value')],
              Input('separator', 'value'),
              )
def scatter_it(analysis, done, display, sep):
    if analysis==1 and done is not None and display==1:
        if sep is None:
            sep='##'
        dft=pd.read_parquet('pre-matrix.gzip')
        df_labels=pd.read_parquet('df_cluster.gzip')
        dft['label']=df_labels
        try:
            dft['number']=dft.name.map(lambda x: look_for_number(x, sep)) 
        except:
            dft.name=dft.name.map(lambda x: str(x)+sep+'200' if sep not in x else x)
            dft['number']=dft.name.map(lambda x: look_for_number(x, sep))
            print("there was an error! Be carreful, every sequence name must contain the separator. A default '200' value was attributed to missing data")
        pc=pd.read_parquet('PC_plot_3D.qzip')
        fig=px.scatter(pd.concat([dft, pc], axis=1),x='PC1', y='PC2', color='file', size='number', 
                       size_max=(pc.PC1.max()-pc.PC1.min())/6, hover_name='name',
                       hover_data=['name', 'label', 'number'])
        fig.add_trace(go.Histogram2dContour(x=pc.PC1, y=pc.PC2, ncontours=80, hoverinfo='skip'))
        fig.update_layout(height=1000)
        fig.update_layout(template='plotly_white')
        return [fig, fig]
    else: 
        return dash.no_update
    """
    [
            html.Div([
            dcc.Graph(figure=fig)
            ])
            ,
            html.Div([
            dcc.Graph(figure=fig)
            ])
            ]
    
    
    [
            html.Div([
            dcc.Graph()
            ])
            ,
            html.Div([
            dcc.Graph()
            ])
            ]"""
    

@app.callback(Output('reads_count', 'children'),
              Input('scatter_contour', 'figure'),
              Input('separator', 'value'),
              )
def count_reads(figure_done, sep):
    if sep is None: 
        sep='##'
    if figure_done is not None:
        dft=pd.read_parquet('pre-matrix.gzip')
        df_labels=pd.read_parquet('df_cluster.gzip')
        dft['label']=df_labels
        try:
            dft['number']=dft.name.map(lambda x: look_for_number(x, sep))
        except:
            dft.name=dft.name.map(lambda x: str(x)+sep+'200' if sep not in x else x)
            dft['number']=dft.name.map(lambda x: look_for_number(x, sep))
        d=dft.groupby(by='file').number.sum()
        return dcc.Graph(
            figure=px.pie(d, values=d.values, names=d.index)
            )

@app.callback(Output('bar_clusters', 'children'),
              Input('scatter_contour', 'figure'),
              Input('order', 'value'),
              Input('upload-data', 'filename'),
              Input('separator', 'value'),
              )
def draw_bar_clusters(figure, order, list_of_names, sep):
    if figure is not None and order is not None:
        if sep is None:
            sep='##'
        #____import datatable
        df_label=pd.read_parquet('df_cluster.gzip')
        df=pd.read_parquet('pre-matrix.gzip')
        df['label']=df_label.label
        #now i get the number of time each sequence were found in each file
        try:
            df['number']=df.name.map(lambda x: look_for_number(x, sep))
        except:
            df.name=df.name.map(lambda x: str(x)+sep+'200' if sep not in x else x)
            df['number']=df.name.map(lambda x: look_for_number(x, sep))
        #___draw bar chart
        fig=px.bar(df, x='numbe', y='label', color='file',template='seaborn')
        return(dcc.Graph(figure=fig))
        
        
        

@app.callback(Output('selected_points_table', 'children'),
              Input('scatter_contour2', 'selectedData'),
              Input('separator', 'value'),
              )
def display_selected(selectedData, sep):
    if sep is None:
        sep='##'
    if selectedData and selectedData['points']:
        selection=[p['hovertext'] for p in selectedData['points']]
        df=pd.read_parquet('pre-matrix.gzip')
        df.name=df.name.map(lambda x: str(x)+sep+'200' if sep not in x else x) #so that changed name can be selected
        dfs=df.loc[[i in selection for i in df.name]]
        
        return dash_table.DataTable(
            data=(dfs.to_dict('records')),
            columns=[{"name": i, "id": i} for i in df.columns],
            page_size=10,
            export_format="csv",
            )


@app.callback(Output('cluster_table', 'children'),
              Input('selecte_cluster', 'value'),
              )
def display_cluster(cluster):
    if cluster is not None:
        df_label=pd.read_parquet('df_cluster.gzip')
        df=pd.read_parquet('pre-matrix.gzip')
        df['label']=df_label.label
        d=df.loc[df.label==cluster]
        return dash_table.DataTable(
            data=(d.to_dict('records')),
            columns=[{"name": i, "id": i} for i in d.columns],
            page_size=10,
            export_format="csv",
            )

@app.callback(Output('cluster_pie', 'children'),
              Input('selecte_cluster', 'value'),
              Input('scatter_contour2', 'selectedData'),
              Input('separator', 'value'),
              )
def cluster_pie(cluster, selectedData, sep):
    if sep is None:
        sep='##'
        print('A separator must be given, otherwise "##" isgiven as a default value and every sequence are counted 200 times')
    #_________Here we plot the pie chart of file contribution to selected cluster. It has prioroty over downstream code
    if cluster is not None:
        df_label=pd.read_parquet('df_cluster.gzip')
        df=pd.read_parquet('pre-matrix.gzip')
        df['label']=df_label.label
        df=df.loc[df.label==cluster]
        #now i get the number of time each sequence were found in each file
        try:
            df['number']=df.name.map(lambda x: look_for_number(x, sep))
        except:
            df.name=df.name.map(lambda x: str(x)+sep+'200' if sep not in x else x)
            df['number']=df.name.map(lambda x: look_for_number(x, sep))
        d=df.groupby(by='file').sum()
        return dcc.Graph(
            figure=px.pie(d, values=d.number, names=d.index)
            )
    #_________Here we plot the pie chart of file contribution to selected point on the scatter plot
    elif selectedData and selectedData['points']:
        selection=[p['hovertext'] for p in selectedData['points']] #get selected
        df=pd.read_parquet('pre-matrix.gzip') #import data
        df_label=pd.read_parquet('df_cluster.gzip') #import cluster attribution
        df['label']=df_label.label #add label
        df.name=df.name.map(lambda x: str(x)+sep+'200' if sep not in x else x) #correct hover that dont' contain the separator
        df=df.loc[[i in selection for i in df.name]] #keep selected points
        df['number']=df.name.map(lambda x: look_for_number(x, sep)) #add numbers so we can calculate the file's contribution
        d=df.groupby(by='file').sum()
        return dcc.Graph(
            figure=px.pie(d, values=d.number, names=d.index)
            )
        
@app.callback(Output('cluster_plot', 'children'),
              Input('selecte_cluster', 'value'),
              Input('separator', 'value'),
              )
def plot_selected_cluster(cluster, sep):
    if cluster is not None:
        if sep is None:
            sep='##'
        dft=pd.read_parquet('pre-matrix.gzip')
        df_labels=pd.read_parquet('df_cluster.gzip')
        dft['label']=df_labels
        dft=dft.loc[dft.label==cluster]
        split=split_df(dft)
        c=get_consensus(split)
        m=matrix_df(split, c)
        pca_plot = PCA(n_components=2)
        a=m.columns.map(lambda x: type(x)==int) #only get columns that are AA positions
        pc = pd.DataFrame(pca_plot.fit_transform(m.loc[:, a]))
        pc.columns=['PC1', 'PC2']
        dft['PC1']=pc.PC1
        dft['PC2']=pc.PC2
        try:
            dft['number']=dft.name.map(lambda x: look_for_number(x, sep)) #attention à changer ça avec quelque chose de bien
        except:
            dft.name=dft.name.map(lambda x: str(x)+sep+'200' if sep not in x else x)
            dft['number']=dft.name.map(lambda x: look_for_number(x, sep))
            print("there was an error! Be carreful, every sequence name must contain the separator. A default '200' value was attributed to missing data")
        fig=px.scatter(dft,x=pc['PC1'], y=pc['PC2'], color='file', size='number', 
                       size_max=(pc.PC1.max()-pc.PC1.min())/6, hover_name='name',
                       hover_data=['name', 'label', 'number'])
        fig.update_layout(height=1000)
        fig.update_layout(plot_bgcolor='#191970')
        #fig.update_layout(template='plotly_dark')
        return dcc.Graph(
            figure=fig
            )


#webbrowser.open_new('http://127.0.0.1:8050/')

    
if __name__ == '__main__':
    app.run_server(debug=True,use_reloader=False)
    

    
