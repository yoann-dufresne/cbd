import os
import matplotlib.pyplot as plt
import sys

def main():
    spercent,stime,cpercent,ctime=dataFromRep(sys.argv[1])
    a,b=average(spercent,stime)
    c,d=average(cpercent,ctime)
    print(a,b)
    plt.plot(a,b,label="successor")
    plt.plot(c,d,label="contain")
    plt.legend(loc='upper center')
    plt.minorticks_on()
    plt.title("for 100000 test with kmer")
    plt.xlabel("percentage")
    plt.ylabel("time(s)")
    plt.savefig(sys.argv[2])

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