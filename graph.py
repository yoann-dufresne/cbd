import os
import matplotlib.pyplot as plt
import sys
import plotly.graph_objects as go


def main():
    spercent,stime,cpercent,ctime=dataFromRep(sys.argv[1])
    fig = go.Figure()
    fig.add_trace(go.Box(
        y=stime,
        x=spercent,
        name='successor',
        marker_color='#3D9970',
        boxpoints=False, # no data points to force whisker to go to the extreme value
    ))
    fig.add_trace(go.Box(
        y=ctime,
        x=cpercent,
        name='contains',
        marker_color='#FF4136',
        boxpoints=False, 
    ))
    fig.update_layout(
    yaxis_title='time',
    boxmode='group' # group together boxes of each percent
    )
    fig.show()
    

def processdata(percent,time):
    ptmp=percent[0]
    epercent=[]
    atime=[]
    i=0
    tmptime=[]
    for p in percent:
        if p==ptmp:
            tmptime.append(time[i])
        else:
            epercent.append(ptmp)
            atime.append(tmptime)
            tmptime=[time[i]]
            ptmp=p
        i+=1
    epercent.append(ptmp)
    atime.append(tmptime)
    return epercent,atime
# take two list(sorted) that contain a percentage and a time and average the time for each percentage present
def average(percent,time):
    ptmp=percent[0]
    sum=0
    nb=0
    avepercent=[]
    avetime=[]
    i=0
    for p in percent:
        if p==ptmp:
            sum+=time[i]
            nb+=1
        else:
            avepercent.append(ptmp)
            avetime.append(sum/nb)
            sum=time[i]
            nb=1
            ptmp=p
            print(p)
        i+=1
    avepercent.append(ptmp)
    avetime.append(sum/nb)

    return avepercent,avetime


#return the usefull data from one file
def dataFromFile(path):
    file = open(path, "r")
    for i in range(4):
        file.readline()
    cs=(file.readline())
    percent=int(file.readline())
    time=float(file.readline())
    return (percent, time,cs)
#return all the sorted data from a directory
def dataFromRep(dirpath):
    files = os.listdir(dirpath)
    spercent=[]
    stime=[]
    cpercent=[]
    ctime=[]
    successor=[]
    contain=[]
    for name in files:
        a=dataFromFile(dirpath+"/"+name)
        if(a[2]=='successor\n'):  
            successor.append((a[0],a[1]))
        else:
            contain.append((a[0],a[1]))
    list.sort(successor)
    list.sort(contain)
    for b,c in successor:
        spercent.append(b)
        stime.append(c)
    for b,c in contain:
        cpercent.append(b)
        ctime.append(c)

    return (spercent, stime,cpercent,ctime)






if __name__ == "__main__":
    main()