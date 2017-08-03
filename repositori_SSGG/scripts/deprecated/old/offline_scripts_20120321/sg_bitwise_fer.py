from operator import itemgetter, attrgetter

def compare(a,b):
    print a[1],b[1]
    return a[1] < b[1]

def convierteFLAG(flag):
    ret = []
    if flag & 0x1: # 1
        ret.append(['Paired',1])

        if flag & 0x2: # 2
            ret.append(['Mapped in a proper pair',2])

        if flag & 0x8: # 8
            ret.append(["Mate unmapped",5])
        else:
            if flag & 0x20: # 32
                ret.append(["Mate Strand = Reverse",6])
            else: 
                ret.append(["Mate Strand = Forward",6])

        if flag & 0x40: # 64
            ret.append(["First read of the pair",7])
            
        if flag & 0x80: # 128
            ret.append(["Second read of the pair",8])
    
    if flag & 0x4: # 4
        ret.append(['Sequence unmapped',3])

    if flag & 0x10: # 16
        ret.append(['Reverse Strand',4])
    else:
        ret.append(['Forward Strand',4])

    if flag & 0x100: # 256
        ret.append(['The alignment is not primary (a read having split hits may have multiple primary alignment records)',9])

    if flag & 0x200: # 512
        ret.append(['The read fails platform/vendor quality checks',10])

    if flag & 0x400: # 1024
        ret.append(['The read is either a PCR duplicate or an optical duplicate',11])
    #print sorted(ret,key=itemgetter(1))

    ret.sort(key=itemgetter(1))
    ret = [v[0] for v in ret]
    return ret

import sys

cad = convierteFLAG(int(sys.argv[1]))
print cad
