import os
import matplotlib.pyplot as plt

def main():
    a,b=dataFromRep("resultpercentseq")
    
    plt.plot(a,b)
    plt.minorticks_on()
    plt.show()


def dataFromFile(path):
    file = open(path, "r")
    for i in range(5):
        file.readline()
    percent=int(file.readline())
    time=float(file.readline())
    return (percent, time)

def dataFromRep(dirpath):
    files = os.listdir(dirpath)
    percent=[]
    time=[]
    tmp=[]
    for name in files:
        a=dataFromFile(dirpath+"/"+name)
        tmp.append(a)
    list.sort(tmp)
    for b,c in tmp:
        percent.append(b)
        time.append(c)
    return [percent, time]






if __name__ == "__main__":
    main()