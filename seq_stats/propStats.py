#run with python3 propStats.py

#For a file containing proportion of sequenced reads successfully aligned to the genome (generate with addPropAligned.sh), calculate the average, min and max proportion across samples

fhand=open('propReadsAligned.txt')
props=list()
for line in fhand:
    line=line.rstrip()
    words=line.split()
    prop=words[3]
    try:
        fprop=float(prop)
    except:
        continue
    else:
        props.append(fprop)
#print(props)
print('Average:',sum(props)/len(props))
print('Min:',max(props))
print('Max:',min(props))
