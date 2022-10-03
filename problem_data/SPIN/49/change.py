import csv
import sys 
data = []
with open ('49_1') as f : 
    for line in f : 
        tt = line.split(' ')
        temp =[]
        s = ""
        if(len(tt)>1) :
            temp.append(int(tt[0]))
            temp[0]+= 1
            temp.append(int(tt[1]))
            temp[1] += 1
            temp.append(tt[2])
            s += str(temp[0])
            s += ' '
            s += str(temp[1])
            s += ' '
            s += temp[2]
        else  :
            s += tt[0] 
        data.append(s)

ww = open ('49_1','w')
for line in data :
   ww.write(line)
ww.close()


