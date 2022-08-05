import os
import matplotlib.pyplot as plt
import sys
import plotly.graph_objects as go


def main():
    divider=testnumber(sys.argv[1])
    spercent,stime,cpercent,ctime=dataFromRep(sys.argv[1],divider)
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
    xaxis_title='pourcentage de sequence connus',
    boxmode='group', # group together boxes of each percent
    title="temps en seconde pour un kmer"
    )
    fig.write_image(sys.argv[1]+".png", width=1000, height=500, scale=2)
    
def testnumber(path):
    files = os.listdir(path)
    file = open(path+"/"+files[0], "r")
    for i in range(2):
        file.readline()
    return int(file.readline())



#return the usefull data from one file
def dataFromFile(path,divider):
    file = open(path, "r")
    for i in range(4):
        file.readline()
    cs=(file.readline())
    percent=int(file.readline())
    time=float(float(file.readline())/divider)
    return (percent, time,cs)
#return all the sorted data from a directory
def dataFromRep(dirpath,divider):
    files = os.listdir(dirpath)
    spercent=[]
    stime=[]
    cpercent=[]
    ctime=[]
    successor=[]
    contain=[]
    for name in files:
        a=dataFromFile(dirpath+"/"+name,divider)
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