#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np

#Extraction of data from perfIsPresent.txt
def getDataFromFile():
    file = open('perfIsPresent.txt', "r")
    title = [] #Will represent each percentages from 0 to 100
    value = [] #Will represent time value of each percentage
    for line in file:
        sp = line.split('\t') #data are separate with a '\t'
        title.append(int(sp[0])) #first column is percentages
        value.append(float(sp[1])) #second is values
    return [title, value]

#Will build graphs depending on data
#We build a normal graph (plot) to have an overall view of the function behaviour
#and a bar graph to see evolution more precisely
def graphsConstructor():
    dataForGraphs = getDataFromFile() #data recovery
    fig, axs = plt.subplots(2) #We will plot 2 graphs in the same time
    fig.suptitle("Evalutation of time performance for isPresent")
    axs[0].plot(dataForGraphs[0], dataForGraphs[1]) #plot graph
    axs[0].set_xlabel('Presence in percentage')
    axs[0].set_ylabel('Time in microseconds')
    axs[1].bar(dataForGraphs[0], dataForGraphs[1])  #bar graph
    axs[1].set_xlabel('Presence in percentage')
    axs[1].set_ylabel('Time in microseconds')
    axs[1].grid(color='#95a5a6', linestyle='--', linewidth=2, axis='y', alpha=0.7)
    plt.show()

def main():
    graphsConstructor()
    
if __name__ == "__main__":
    main()
