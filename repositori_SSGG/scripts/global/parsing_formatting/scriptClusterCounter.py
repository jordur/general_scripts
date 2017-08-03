#!/share/apps/local/python/bin/python

import sys

inputFilename = 'out.extendedAssembledFrags.fasta.clstr'
outputFilename = 'clusterCounter'

try:
    inputFilename = sys.argv[1]
    outputFilename = sys.argv[2]
except:
    print ("The call for this script is: python scriptClusterCounter.py inputFilename outputFilename")
    exit(1)


try:
    #f = open('out.extendedAssembledFrags.fasta.clstr', 'r')
    #g = open('clusterCounter', 'w')
    f = open(inputFilename, 'r')
    g = open(outputFilename, 'w')
except:
    print ("Error opening files.")
    exit(1)

actualClusterNumber = 0
actualClusterName = "vacio"
firstLine = True

for line in f:
    #print line,

    firstChar = line[0]
    comparison = firstChar is '>'

	#if firstChar in '>':
    if comparison:
        #It is a new cluster
        if firstLine:
            g.write('ClusterName\tCounter\n')
            firstLine = False
        else:
            g.write(actualClusterName+'\t'+str(actualClusterNumber)+'\n')

        actualClusterNumber = 0
        actualClusterName  = line[1:-1]     #Es un substring quitando el primer caracter y el ultimo pk es un intro
        actualClusterName  = actualClusterName.replace(' ','')
        #actualClusterName = line[1:-2]		#Nos quita los dos ultimos caracteres
        #actualClusterName = line[::-1]     #Nos da la linea invertida!
    else:
        actualClusterNumber += 1

f.close()
g.close()

