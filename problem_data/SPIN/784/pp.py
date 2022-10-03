

s = "cp /home/andy3466/FastEfficientP3/problem_files/IsingSpinGlass_pm_784_"  

for i in range(101): 
    if i == 0 : 
        continue 
    temp = s 
    temp += str(i) 
    temp += ".txt 784_"
    temp += str(i)
    print temp 
    print "python change.py %d" %(i)


