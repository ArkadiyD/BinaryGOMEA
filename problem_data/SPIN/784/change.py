import csv
import sys 
data = []
datav = []
s = "784_" 
s += sys.argv[1]
with open (s) as f : 
    for line in f :
        data.append(line)
tx = int(data[1])
tx /= 2 
tk = str(tx)
tk +='\n'
datav.append(tk)
tx = data[0].split(' ')
ttt = float(tx[0])
ttt /= 784
ttt = round(ttt,5)
ttt += 0.00002
tk = str(ttt)
tk += '\n'
datav.append(tk)

for i in range(len(data)):
    if i > 1 : 
       temp = data[i].split(' ')
       ttt =""
       tt = int(temp[0])
       tt += 1
       ttt += str(tt)
       ttt += " "
       tt = int(temp[1])
       tt +=1
       ttt += str(tt)
       ttt += " "
       ttt += temp[2]
       datav.append(ttt)

ww = open (s,'w')
for line in datav :
   ww.write(line)
ww.close()


