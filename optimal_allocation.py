#!/usr/bin/env python
# coding: utf-8

#Citation Request:

#Please cite the article: Shima Azizi, Cem Deniz Caglar Bozkir, Andrew C. Trapp, O. Erhun Kundakcioglu, Ali Kaan Kurbanzade, "Aid Allocation for Camp-Based and Urban Refugees with Uncertain Demand and Replenishments", 2020.

#Email: sazizi@wpi.edu, cem.bozkir@ozu.edu.tr, atrapp@wpi.edu, erhun.kundakcioglu@ozyegin.edu.tr, kaan.kurbanzade@ozu.edu.tr

import pandas as pd
import numpy as np
import math
import sys
import os
from pulp import *
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import time
from operator import add 

Data = pd.read_csv('rates.csv', encoding='utf-8',delimiter=',')
df2 = pd.read_csv('parameters.csv', encoding='utf-8', delimiter=',')
h = df2.iat[0,0]
mu = 1 / df2.iat[0,1]
deltar = df2.iat[0,2]
deltad = df2.iat[0,3]
alpha = df2.iat[0,4]

if alpha>mu:
    #print("Deprivation Rate cannot exceed yearly cycle rate.")
    sys.exit("Deprivation Rate cannot exceed yearly cycle rate.")

if deltar>(deltad*(alpha/(mu-alpha))):
    sys.exit("The ratio between Referral Cost Per External Demand and Deprivation Cost Per Time Unit Per Internal Demand does not hold.")
    

Data=Data.rename(columns={"Camp Name": "i", "Internal Rate of Arrival (/yr)": "lambda_C","External Rate of Arrival (/yr)": "lambda_S"})

Data['omega'] = math.log(deltar/(deltad*(alpha/(mu-alpha))))/np.log(Data['lambda_C']/(Data['lambda_C']+mu))
Data['gamma1_1'] = (((Data['lambda_C']+Data['lambda_S']+mu)/(Data['lambda_C']+Data['lambda_S']))**Data['omega'])*(Data['lambda_S']*(h/(mu**2)+deltar/mu)+(Data['lambda_C']*((Data['lambda_C']/(Data['lambda_C']+mu))**Data['omega'])*(h/(mu*mu)+deltad*(alpha/(mu-alpha)))))
Data['gamma1_2'] = (Data['lambda_C']+Data['lambda_S'])/(Data['lambda_C']+Data['lambda_S']+mu)
Data['gamma1_3'] = h/mu
Data['gamma1_4'] = (h/(mu**2))*(-Data['lambda_C']-Data['lambda_S'])
Data['gamma2_1'] = deltad*Data['lambda_C']*(alpha/(mu-alpha))+h*Data['lambda_C']/(mu*mu)
Data['gamma2_2'] = (Data['lambda_C'])/(Data['lambda_C']+mu)
Data['gamma2_3'] = h/mu
Data['gamma2_4'] = Data['lambda_S']*deltar/mu-h*Data['lambda_C']/(mu*mu)


Data.insert(3, 'h', h)
Data.insert(4, 'Mu', mu)
Data.insert(5, 'delta_R', deltar)
Data.insert(6, 'delta_D', deltad)
Data.insert(7, 'alpha', alpha)


def get_files(dirname, reverse=False):
    instances = []      
    for file_name in os.listdir(path):
        if file_name.endswith('.csv'):
            instances.append(os.path.join(path,file_name))            
    SortedInstances=sorted(instances)    
    return SortedInstances 

def ExportResult(ListofResults,N):
        cols = ['Problem Code',
                'Supply',
                'Sum over Thresholds',
                'Overall Demands-Inside',
                'Overall Demands-Outside',
                'Overall Demands',
                'Step',
                'Number of breakpoints']
        
        for i in range(N):
            cols.append("Camp "+str(i+1))
            
        cols = cols + ['Sum over X Vector',
                       'Model Runtime',
                       'Objective Value'                              
                         ]  
        df = pd.DataFrame(ListofResults,columns = cols)
        path = os.path.join(os.getcwd(),'Results.csv')
        df.to_csv(path,index=False)

def tickformat(x): #function for customize axis settings
    if int(x) == float(x):
        return str(int(x))
    else:
        return str(x)      
    
def PWL(Data,problem_Code):
    N = len(list(Data['i']))     
    I = [0] * N  # List of initial inventory for each camp
    threshold = list(Data['omega'])     
    Demand_Inside = list(Data['lambda_C'])   
    Demand_Outside = list(Data['lambda_S'])   
    gamma1_1 = list(Data['gamma1_1'])    
    gamma1_2 = list(Data['gamma1_2'])        
    gamma1_4 = list(Data['gamma1_4'])    
    gamma2_1 = list(Data['gamma2_1'])    
    gamma2_2 = list(Data['gamma2_2'])        
    gamma2_4 = list(Data['gamma2_4']) 
    h= list(Data['h'])[0] 
    Mu = list(Data['Mu'])[0]
    Sum_over_Threshold = sum(threshold)
    Total_Demand_Inside = sum(Demand_Inside)
    Total_Demand_Outside = sum(Demand_Outside)
    Total_Demand = Total_Demand_Inside + Total_Demand_Outside
    SList = list(map(int,[Sum_over_Threshold/10 , Sum_over_Threshold/5 , Sum_over_Threshold/3,  Sum_over_Threshold/2]))
    SList.extend(list(map(int, np.linspace(Sum_over_Threshold,max(Sum_over_Threshold*3, Total_Demand), 5)))) 


    for S in SList:   
        UB_x = [S] * N
        LB_x = [0] * N 
        if S <= 50000:
            Step =2               
        elif 50000 < S <= 500000:
            Step =10
        elif 500000 < S :
            Step=50
        
        Number_of_Breakpoints = S//Step + 1 + 1
        D_Below=[]  #List of number of breakepoints for each camp for below threshold
        D_Above=[]  #List of number of breakepoints for each camp for above threshold
        for i in range(N):
            if S > threshold[i]:
                D_Below.append((threshold[i]//Step)+1)
                D_Above.append(Number_of_Breakpoints - ((threshold[i]//Step)+1) - 1)            
            elif S <= threshold[i]:
                D_Below.append((S//Step)+1)
                D_Above.append(0)            
            
        # Create Breakpoints:
        B = {}       # b_i^j
        T_B = {}     # T(b_i^j)         
        for i in range(N):  
            First_B_After_Threshold = S - (Number_of_Breakpoints-( D_Below[i] +2)) * Step
            for j in range(Number_of_Breakpoints):
                if j < D_Below[i]: 
                    B[i,j] = LB_x[i] + Step * (j)  
                elif j == D_Below[i]:
                    B[i,j] = threshold[i]     
                else:
                    B[i,j] = First_B_After_Threshold + Step * (j- D_Below[i] - 1)  
                                    
        for i in range(N):
            for j in range(Number_of_Breakpoints):                 
                    if B[i,j] > threshold[i]:
                        T_B[i,j] = (gamma1_1[i] * (gamma1_2[i] ** B[i,j])) +  gamma1_4[i]          
                    elif B[i,j] <= threshold[i]: 
                        T_B[i,j] = (gamma2_1[i] * (gamma2_2[i] ** B[i,j]))  + gamma2_4[i]   
                    else:
                        return "Error!"
                                                 
        # compute slope of the j-th segment:
        m = {}
        c = {}                       
        for i in range(N):
            for j in range(1,Number_of_Breakpoints):
                    m[i,j] = (T_B[i,j] - T_B[i,j-1]) / (B[i,j] - B[i,j-1])                       
                    c[i,j] =  T_B[i,j] - (m[i,j] * B[i,j])                                             

        m1 = LpProblem('PWL-MCH', LpMinimize)
        # Create variables
        x = {}  # initial inventory for camp i
        xi = {} # x_i^j = x[i][j] if x[i][j] lies in the j-th segment for camp i; and  0 otherwise
        g = {}  # g_i^j # g_i^j is a binary variable, g[i,j] = 1 if x[i,j] lies in j segment; and 0 otherwise. 
        nu = {}             
        for i in range(N):
            x[i] = LpVariable(lowBound = 0, upBound = S, cat = 'Continuous', name='x(%d)'%(i) )
            nu[i] = LpVariable(lowBound = -1e20, upBound = 1e20, cat = 'Continuous', name='nu(%d)'%(i) )
            for j in range(1,Number_of_Breakpoints):
                xi[i,j] = LpVariable(lowBound = 0, upBound = S, cat = 'Continuous', name='xi(%d,%d)'%(i,j) )
                g[i,j] = LpVariable(lowBound = 0.0, upBound = 1.0, cat = 'Binary',name='g(%d,%d)'%(i,j))                

        # Constraints:    
        for i in range(N):       
            m1 += lpSum([xi[i,j] for j in range(1,Number_of_Breakpoints)]) == x[i],'c1(%d)'%(i)

            m1 += lpSum([(m[i,j]*xi[i,j]) + (c[i,j]*g[i,j]) for j in range(1,Number_of_Breakpoints)]) == nu[i],'c2(%d)'%(i)

            m1 += lpSum([g[i,j] for j in range(1,Number_of_Breakpoints)]) == 1,'c3(%d)'%(i)
            for j in range(1,Number_of_Breakpoints):
                m1 += B[i,j-1]*g[i,j] <= xi[i,j], 'c4(%d,%d)'%(i,j)
                m1 += xi[i,j] <= B[i,j]*g[i,j], 'c5(%d,%d)'%(i,j)           
        m1 += lpSum([(x[i] - I[i]) for i in range(N)]) <= S, 'c6' 
           
        #Define objective function
        m1 += lpSum([nu[i] for i in range(N)]) + (h/Mu)*(S+lpSum([I[i] for i in range(N)]))
        start_time = time.time() 
        m1.solve(PULP_CBC_CMD())
        Execution_time =  time.time() - start_time
        l = [problem_Code[:-4],
                             S,
                             int(Sum_over_Threshold),
                             int(Total_Demand_Inside),
                             int(Total_Demand_Outside),
                             int(Total_Demand),
                             Step,
                             Number_of_Breakpoints]
        X_Vector=[]
        for i in range(N):
                l.append(x[i].varValue)
                X_Vector.append(x[i].varValue)             
        l = l + [sum(X_Vector),
                 Execution_time,
                 value(m1.objective)]
        ListofResults.append(l)

    ExportResult(ListofResults,N)


def Pie_Demands_Seperate(Data,problem_Code): 
 
        N = len(list(Data['i'])) # Number of camps
        Demand_Inside = list(Data['lambda_C'])    
        Demand_Outside = list(Data['lambda_S'])        
        
        Total_Demand_List = list(map(add, Demand_Inside, Demand_Outside))   #list of total demand for each camp

        idx_Total_Demand_List = sorted(range(len(Total_Demand_List)), key=lambda k: Total_Demand_List[k], reverse=True) #index of sorted (largest to smallest demand)
        
        r = 5*0.3
        fig, ax = plt.subplots(nrows=N, ncols=1)

        row=0
        for i in idx_Total_Demand_List:

            ax[row].pie([Demand_Inside[i],Demand_Outside[i]],radius = r, shadow=False,colors=["lightcoral","lightskyblue"],wedgeprops={"edgecolor":"white",'linewidth': 1, 'antialiased': True})
            ax[row].set_title("Camp "+str(i+1),fontsize=16)
            ax[row].title.set_position([2, 0.3])

            r -= 0.08
            row += 1
 
        labels=["Internal Demand Rate","External Demand Rate"]
        ax[0].legend(labels,bbox_to_anchor=(1, 1), loc=8,fontsize=16)
        fig = plt.gcf()
        fig.set_size_inches(10, 8)
        path_plot = os.path.join(os.getcwd(),'Plots_'+problem_Code+"_Demands")
        fig.savefig(path_plot, dpi=300,bbox_inches = 'tight') 

def Rate_Above_Threshold(Result,Data,problem_Code):
    N = len(list(Data['i'])) # Number of camps
    Demand_Inside = list(Data['lambda_C'])    
    Demand_Outside = list(Data['lambda_S'])  
    Total_Demand = [Demand_Inside[i]+Demand_Outside[i] for i in range(N)] 
    Supply=list(Result["Supply"])
    Supply_Max= max(Result["Supply"])
    Index_Max_Supply = Supply.index(max(Supply))
    
    Label=[]
    for i in range(N):
        Label.append("Camp "+str(i+1))
        
    #Check if it is above supply:
    if Result["Supply"][Index_Max_Supply-2] > Result["Sum over Thresholds"][0]:
        Supply_Above_Threshold = Result["Supply"][Index_Max_Supply-2]    
    else:
        return "Error"


    idx_Total_Demand = sorted(range(len(Total_Demand)), key=lambda k: Total_Demand[k], reverse=True) #index of sorted (largest to smallest demand)
    #Sort all demands based on total demands:
    Demand_Inside = [Demand_Inside[i] for i in idx_Total_Demand]
    Demand_Outside = [Demand_Outside[i] for i in idx_Total_Demand]
    Total_Demand = [Total_Demand[i] for i in idx_Total_Demand]
    Label = [Label[i] for i in idx_Total_Demand]
    #Calcute Rates:
    Denominator = (Supply_Max)-(Supply_Above_Threshold) 
    List_Rate_Above=[]
    for i in idx_Total_Demand:
       Rate= (Result["Camp "+str(i+1)][Index_Max_Supply] - Result["Camp "+str(i+1)][Index_Max_Supply-2])/Denominator
       List_Rate_Above.append(Rate)                       
#Right axis:
    fig, ax1 = plt.subplots()

    I = ax1.plot(Label, Demand_Inside, color="black", linewidth=2, alpha=0.9,linestyle='--', label="Internal Demand Rate")
    O = ax1.plot(Label,Demand_Outside, color="black", linewidth=2, alpha=0.9,linestyle=':', label="External Demand Rate")
    T = ax1.plot(Label,Total_Demand, linewidth=6, alpha=0.9,linestyle='-', label="Overall Demand Rate")
    ax1.tick_params(axis='y')
    plt.yticks(fontsize=20)    
    plt.xticks(fontsize=15,rotation=45) 
#Left Axis    
    ax2 = ax1.twinx()   
    R=ax2.scatter(Label,List_Rate_Above, color="red", marker='s',label = "Fraction of Allocated Aid")    
    ax2.tick_params(axis='y')  
    fig.tight_layout()     
    plt.setp(ax1.spines.values(), linewidth=0.1)
    plt.setp(ax2.spines.values(), linewidth=0.1)        
    #format x and y axis to be in thousand
    fmt = FuncFormatter(lambda x, pos: tickformat(x / 1e3))
    ax1.yaxis.set_major_formatter(fmt)
    ax2.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y)))     
    #Add legend:
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper right',fontsize=20)

    plt.yticks(fontsize=20) 
    plt.title("Optimal Allocation of Aid as the Fraction of Supply Allocated \n to Each Camp - Above Transition Zone.",fontsize=20)          
    fig = plt.gcf()
    fig.set_size_inches(15, 10)#for 10 Camp
    path_plot = os.path.join(os.getcwd(),'Plots_'+problem_Code+"_Rate_Above")
    fig.savefig(path_plot, dpi=300,bbox_inches = 'tight')
    #plt.show()


def Rate_Below_Threshold(Result,Data,problem_Code):   
    N = len(list(Data['i'])) # Number of camps
    Demand_Inside = list(Data['lambda_C'])    
    Demand_Outside = list(Data['lambda_S'])  
    Total_Demand = [Demand_Inside[i]+Demand_Outside[i] for i in range(N)] 
    Supply=list(Result["Supply"])
    Label = []
    for i in range(N):
        Label.append("Camp "+str(i+1))
    
    Supply_Min= min(Result["Supply"])
    Index_Min_Supply = Supply.index(min(Supply))

    #Check if it is below supply:
    if Result["Supply"][Index_Min_Supply+1] < Result["Sum over Thresholds"][0]:
        Supply_Below_Threshold = Result["Supply"][Index_Min_Supply+1]    
    else:
        return "Error"


    idx_Demand_Inside = sorted(range(len(Demand_Inside)), key=lambda k: Demand_Inside[k], reverse=True) #index of sorted (largest to smallest demand)
    #Sort all demands based on total demands:
    Demand_Inside = [Demand_Inside[i] for i in idx_Demand_Inside]
    Demand_Outside = [Demand_Outside[i] for i in idx_Demand_Inside]
    Total_Demand = [Total_Demand[i] for i in idx_Demand_Inside]
    Label = [Label[i] for i in idx_Demand_Inside]
    #Calcute Rates:
    Denominator = (Supply_Below_Threshold)  - (Supply_Min)
    List_Rate_Below=[]
    for i in idx_Demand_Inside:
       Rate= (Result["Camp "+str(i+1)][Index_Min_Supply+1] - Result["Camp "+str(i+1)][Index_Min_Supply])/Denominator
       List_Rate_Below.append(Rate)
        
#Right axis:
    fig, ax1 = plt.subplots()
#    ax1.set_xlabel('Camp Name',fontsize=20)
#    ax1.set_ylabel('Demand Rate (in thousands)',fontsize=18)
    I = ax1.plot(Label, Demand_Inside, color="black", linewidth=6, alpha=0.9,linestyle='--', label="Internal Demand Rate")
    O = ax1.plot(Label,Demand_Outside, color="black", linewidth=2, alpha=0.9,linestyle=':', label="External Demand Rate")
    T = ax1.plot(Label,Total_Demand, linewidth=2, alpha=0.9,linestyle='-', label="Overall Demand Rate")
    ax1.tick_params(axis='y')
    plt.yticks(fontsize=20)    
    plt.xticks(fontsize=15,rotation=45) 

#Left Axis    
    ax2 = ax1.twinx()  
    
    R=ax2.scatter(Label,List_Rate_Below, color="red", marker='s',label = "Fraction of Allocated Aid")
    
    ax2.tick_params(axis='y')
    
    fig.tight_layout()  
    
    
    plt.setp(ax1.spines.values(), linewidth=0.1)
    plt.setp(ax2.spines.values(), linewidth=0.1)    
    
    #format x and y axis to be in thousand
    fmt = FuncFormatter(lambda x, pos: tickformat(x / 1e3))
    ax1.yaxis.set_major_formatter(fmt)
    ax2.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y))) 
    
    #Add legend:
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper right',fontsize=20)


    plt.yticks(fontsize=20)        
    plt.title("Optimal Allocation of Aid as the Fraction of Supply Allocated \n to Each Camp - Below Transition Zone.",fontsize=20)      
    fig = plt.gcf()
    fig.set_size_inches(15, 10)  
    path_plot = os.path.join(os.getcwd(),'Plots_'+problem_Code+"_Rate_Below")
    fig.savefig(path_plot, dpi=300,bbox_inches = 'tight')
    #plt.show()


def Analyze_X_Spagetti(Result,Data,problem_Code):
    N = len(list(Data['i'])) # Number of camps
    threshold = list(Data['omega'])
    if N<=10:
        Color= ["r","sienna","darkgreen","darkcyan","royalblue","mediumvioletred","goldenrod","dimgrey","chocolate","deeppink"]
    else:
        Color=["steelblue"] * N
    fig,axes=plt.subplots()                    
#Add Zone of Ambiguity:         
    List_Supply_Hit_Threshold = []
    for i in range(N):   
        for j in range(len(Result['Supply'])):          
                if Result["Camp "+str(i+1)][j] > Data['omega'][i]: 
                    List_Supply_Hit_Threshold.append(Result['Supply'][j-1])
                    List_Supply_Hit_Threshold.append(Result['Supply'][j])
#                    List_Supply_Hit_Threshold.append(Result['Supply'][j+1])
                    break
    plt.axvspan(min(List_Supply_Hit_Threshold), max(List_Supply_Hit_Threshold), facecolor = "darkgrey",edgecolor='dimgrey', alpha=0.2,lw=0.5) #Grey:112        

    num=0
    for i in range(N):        

# Add pointer for threshold
        for j in range(len(Result['Supply'])):          
                if Result["Camp "+str(i+1)][j] > Data['omega'][i]: 
                    m = (Result["Camp "+str(i+1)][j] - Result["Camp "+str(i+1)][j-1])/(Result["Supply"][j]-Result["Supply"][j-1])
                    x = (Data['omega'][i]-Result["Camp "+str(i+1)][j] + m * Result["Supply"][j])/m

                    plt.scatter(x, Data['omega'][i], marker='o', s=30, edgecolor=Color[i], c="none", alpha=0.9,linewidth=1)
                    break  
     
        Below=[]
        Above=[]
        for j in range(len(Result['Supply'])):
            if Result["Camp "+str(i+1)][j] <= Data['omega'][i]:   
                Below.append(j)
            else:
                Above.append(j)
        
        plt.plot(Result['Supply'][:max(Below)+1], Result["Camp "+str(i+1)][:max(Below)+1], marker='x', markersize=2, color=Color[num], linewidth=1.5, alpha=0.9)
        plt.plot(Result['Supply'][max(Below):min(Above)+1], Result["Camp "+str(i+1)][max(Below):min(Above)+1], color=Color[num], linewidth=1.5, alpha=0.9)
        plt.plot(Result['Supply'][min(Above):], Result["Camp "+str(i+1)][min(Above):], marker='x', markersize=5, color=Color[num], linewidth=1.5, alpha=0.9, label="Camp "+str(i+1))

        plt.text(max(Result['Supply'])*1.01,max(Result["Camp "+str(i+1)])*0.99,"Camp "+str(i+1),fontsize=14)
        num+=1  
    axes.set_xbound(upper=max(Result['Supply'])*1.1)            
    plt.setp(axes.spines.values(), linewidth=0.1)
    plt.yticks(fontsize=20)    
    plt.xticks(fontsize=20)    
    plt.title("Optimal Allocation Versus Various Levels of Total Supply.",fontsize=20)   
    fmt = FuncFormatter(lambda x, pos: tickformat(x / 1e3))
    axes.xaxis.set_major_formatter(fmt)
    axes.yaxis.set_major_formatter(fmt)   
    plt.xlabel("Total Supply (in thousands)",fontsize=20)
    plt.ylabel("Allocation (in thousands)",fontsize=20)
    fig = plt.gcf()
    fig.set_size_inches(15, 10)
    path_plot = os.path.join(os.getcwd(),'Plots_'+problem_Code+"_Allocation")
    fig.savefig(path_plot, dpi=300,bbox_inches = 'tight')


def Analyze_Optimal_Obj(Result,problem_Code):     

    fig, axes = plt.subplots()
    
    axes.plot(Result["Supply"], Result["Objective Value"], marker='x', markersize=4, linewidth=1.5, alpha=0.9)        
 
    plt.yticks(fontsize=20)    
    plt.xticks(fontsize=20)
    #format x and y axis to be in thousand
    fmt = FuncFormatter(lambda x, pos: tickformat(x / 1e3))
    axes.xaxis.set_major_formatter(fmt)
    axes.yaxis.set_major_formatter(fmt)
    plt.setp(axes.spines.values(), linewidth=0.1) 
    plt.title("Objective Function Versus Total Supply.",fontsize=20)   
    axes.set_ylabel("Expected Total Cost (in thousands)", fontsize=20)
    axes.set_xlabel("Total Supply (in thousands)", fontsize=20)
    fig.set_size_inches(12, 8)
    path_plot = os.path.join(os.getcwd(),'Plots_'+problem_Code+'_Total Expected Cost')
    fig.savefig(path_plot, dpi=300,bbox_inches = 'tight')   
    

ListofResults=[]

problem_Code='single_instance'

PWL(Data,problem_Code)
Pie_Demands_Seperate(Data,problem_Code)
path = os.path.join(os.getcwd())
SortedResults = get_files(path , reverse=False)
Result = pd.read_csv('Results.csv', delimiter=',' )  

Analyze_X_Spagetti(Result,Data,problem_Code)
Rate_Above_Threshold(Result,Data,problem_Code)
Rate_Below_Threshold(Result,Data,problem_Code)
Analyze_Optimal_Obj(Result,problem_Code)
