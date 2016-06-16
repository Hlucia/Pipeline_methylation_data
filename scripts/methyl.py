import numpy
import pandas 
import matplotlib
matplotlib.style.use('ggplot')
import os
def getplot (onlymethyl,nbecotype,min_size,max_size,window) :
    window_list=numpy.arange(min_size,max_size+window,window)
    columnss=numpy.arange(0,nbecotype)
    nbofmethyl=pandas.DataFrame(index=window_list[:-1], columns=columnss)
    for n in range (0, nbecotype):
        line=[]
        for i in range (1, int((max_size+window-min_size)/window)) : 
            line.append(onlymethyl[n].methylation_call[(onlymethyl[n].pos>window_list[i-1])&(onlymethyl[n].pos<=window_list[i])].sum())
        nbofmethyl[n]=line
    nbofmethyl.plot()
    #matplotlib.legend(loc='center right', bbox_to_anchor=(1.0, 0.5)) This line isn't working

def plotclass (CGs,CHGs,CHHs,numecotype,min_size,max_size,window) :
    window_list=numpy.arange(min_size,max_size+window,window)
    columnss=['CGs','CHGs','CHHs']
    nbofmethylclass=pandas.DataFrame(index=window_list[:-1], columns=columnss)
    line1=[]
    line2=[]
    line3=[]
    for i in range (1, int((max_size+window-min_size)/window)) : 
        line1.append(CGs[numecotype].methylation_call[(CGs[numecotype].pos > window_list[i-1]) & (CGs[numecotype].pos <= window_list[i])].sum())
        line2.append(CHGs[numecotype].methylation_call[(CHGs[numecotype].pos > window_list[i-1]) & (CHGs[numecotype].pos <= window_list[i])].sum())
        line3.append(CHHs[numecotype].methylation_call[(CHHs[numecotype].pos > window_list[i-1]) & (CHHs[numecotype].pos <= window_list[i])].sum())
    nbofmethylclass['CGs']=line1
    nbofmethylclass['CHGs']=line2
    nbofmethylclass['CHHs']=line3        
    nbofmethylclass.plot()    

