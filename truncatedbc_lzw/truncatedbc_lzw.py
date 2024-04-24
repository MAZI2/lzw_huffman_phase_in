import math
import os
import subprocess

dictLen = 10923 
# phase-in codes threshold
thr = 5461

"""
dictLen = 5999
thr = 2193
dictLen = 2731 
thr = 1365
"""

# prefill the LZW dictionary
def loadDict(typ, f, dec):
    preDict=[]
    header=[]
    """
    if typ == "t":
        preDict=[
           '.',
           '.TH',
           '.SH',
           '.B',
           '.BR',
           '.I',
           '.P',
           '.SS',
           '.nf',
           '.fi',
           '.TP',
           '\-\-',
           '\-',
           '\fI',
           '\fR',
           '\fB',
           '\fP'
           ]
    el
    """
    if typ == "E":
        header=[
            0x7f,
            0x45,
            0x4c,
            0x46,
            ]
    elif typ == "A":
        preDict=[
           'th',
           'he',
           'in',
           'er',
           'an',
           're',
           'nd',
           'on',
           'en',
           'at',
           'ou',
           'ed',
           'ha',
           'to',
           'or',
           'it',
           'is',
           'hi',
           'es',
           'ng',
           'the',
           'and',
           'ing',
           'her',
           'hat',
           'his',
           'tha',
           'ere',
           'for',
           'ent',
           'ion',
           'ter',
           'was',
           'you',
           'ith',
           'ver',
           'all',
           'wit',
           'thi',
           'tio'
           ]
    elif typ == "L":
        preDict=[
            '((',
            'S1',
            'VP1',
            'VP2',
            'NP',
            'PP'
        ]
    elif typ == "C":
        preDict=[
            '#include',
            'int',
            'static',
            'extern',
            'char',
            'void',
            'unsigned',
            'FILE',
            'free',
            'return',
            'while',
            'if (',
            '==',
            '->',
            '&&',
            '||',
            'break',
            'else',
            'else if',
            'switch',
            'case',
            'return 1;',
            'return 0;'
        ]

    
    slovar=dict()
    if dec==0:
        for i in range(255):
            slovar[chr(i)] = i;

        for i in range(len(preDict)):
            slovar[preDict[i]] = i + 255

        slen = len(slovar)
        if typ == "d" or typ == "E":
            for i in range(255):
                slovar[chr(i+255)] = i + slen;
    else:
        for i in range(255):
            slovar[i] = chr(i);

        for i in range(len(preDict)):
            slovar[i + 255] = preDict[i]

        slen = len(slovar)
        if typ == "d" or typ == "E":
            for i in range(255):
                slovar[i + slen] = chr(i+255)

    return slovar, header

# return a phase-in code bit string of integer 
def toBitString(value, totalValues) :
    length = (totalValues-1).bit_length()
    nextPow2 = (1 << length)
    treshold = nextPow2 - totalValues
    returnString = ""
    if value < treshold:
        length -= 1
    else:
        value += treshold

    returnString = "{0:b}".format(value)
    while len(returnString) < length:
        returnString = "0" + returnString;

    return returnString

# read the whole file and decode phase-in codes
def fromBitString(file, totalValues):
    izhod = []

    buf=""
    thrLen=math.floor(math.log2(totalValues))
    while True:
        temp = buf
        buf+=''.join(format(byte, '08b') for byte in file.read(1))
        if temp == buf and len(buf) > 0:
            num = int(buf[:thrLen], 2)
            if num >= thr:
                b = int(buf[thrLen], 2)
                num = thr+((num-thr)<<1) + b
                buf=buf[thrLen+1:]
            else:
                buf=buf[thrLen:]

            izhod.append(num)

            break

        if len(buf) > thrLen:
            num = int(buf[:thrLen], 2)
            if num >= thr:
                b = int(buf[thrLen], 2)
                num = thr+((num-thr)<<1) + b
                buf=buf[thrLen+1:]
            else:
                buf=buf[thrLen:]
            izhod.append(num)

    return izhod

# concatenates input strings and writes bytes to output file
def writeToFile(string, file, buf):
    if len(string) == 0:
        for i in range(8-len(buf)):
            buf+="0"
    buf+=string

    while len(buf)>=8:
        v = int(buf[0:8], 2)
        file.write(v.to_bytes(1, 'big'))
        buf = buf[8:]

    return buf

# LZW decoding
def dekodiranje(vhod: list, n, typ, out) -> list:
    vhod = fromBitString(vhod, n)

    slovar, header = loadDict(typ, vhod, 1)
    index = len(slovar)

    izhod = []

    N = slovar[vhod[0]];
    izhod.append(N)
    K = N
    for kB in vhod[1:]:
        k = kB
        if k in slovar:
            N = slovar[k]
        else:
            N = K + K[0]

        izhod.append(N)
        if len(slovar) + 1 < dictLen:
            slovar[index] = (K + N[0])
            index += 1

        K = N

    izhod = [i for ele in izhod for i in ele]


    for c in izhod[:-1]:
        out.write(ord(c).to_bytes(1, byteorder='big'))

# LZW encoding
def kodiranje(vhod: list, n, typ, file) -> list:
    slovar, header = loadDict(typ, vhod, 0)

    index = len(slovar)
    #if len(header) != 0:
        
    buf = ""
    N = '' 
    while b:=vhod.read(1):
        z = chr(int.from_bytes(b, byteorder='big'))
        Nz = N + z 
        if Nz in slovar:
            N = Nz
        else:
            buf = writeToFile(toBitString(slovar[N], n), file, buf)

            if len(slovar) + 1 < dictLen:
                slovar[Nz] = index 
                index += 1;
            N = z

    buf = writeToFile(toBitString(slovar[N], n), file, buf)
    writeToFile("", file, buf)

# driver function (input file, output file, mode [encode-0 / decode-1])
# compresses / decompresses the file and returns compression ratio
def naloga2_tekma(dat_vhod: str, dat_izhod: str, nacin: int) -> float:
    inp = open(dat_vhod, "rb")
    outf = open(dat_izhod, "wb")

    # unused file type detection for dictionary prefill
    """
    sub_arg = ""
    if nacin == 0:
        sub_arg = dat_vhod
    else:
        sub_arg = dat_izhod
    p = subprocess.Popen('file -b ' + sub_arg, shell=True, 
                         stdout=subprocess.PIPE, 
                         stderr=subprocess.PIPE,
                         universal_newlines=True
                         )
    out, err = p.communicate()
    typ = str(out[0])
    """
    typ = "d"

    if nacin == 0:
        kodiranje(inp, dictLen, typ, outf)
        outf.close()
        orgs = os.path.getsize(dat_vhod)
        cprs = os.path.getsize(dat_izhod)

        R = orgs / cprs 
    elif nacin == 1:
        dekodiranje(inp, dictLen, typ, outf)
        outf.close()
        cprs = os.path.getsize(dat_vhod)
        orgs = os.path.getsize(dat_izhod)

        R = orgs / cprs 

    return R

# compress 
files = os.listdir("testfiles")

avg=0
count=0
for path in files:
    abs_path = "testfiles/" + path 


    print(path)
    R = naloga2_tekma(abs_path, os.path.splitext(path)[0] + ".lzw", 0)
    R = naloga2_tekma(os.path.splitext(path)[0] + ".lzw", "decoded_" + path, 1)
    print(R)

    """
    original_size = os.path.getsize(abs_path)
    compressed_size = os.path.getsize(os.path.splitext(os.path.basename(abs_path))[0] + ".lzw") 
    decompressed_size = os.path.getsize("decoded_"+os.path.basename(abs_path))
    print(f"Original file size: {original_size} bytes")
    print(f"Decompressed: {decompressed_size/original_size}")
    print(f"Compressed file size: {compressed_size} bytes")
    print(f"Ratio: {compressed_size/original_size}")
    """

    count+=1
    avg+=R

print("Average: ", avg/count)

"""
for k in range(5, 15):
    hi=int(math.pow(2, k))
    lo=int(math.pow(2, k-1))

    maxn=0
    maxr=hi
    maxnum=0
    for n in range(lo, hi):
        num=0
        for i in range(n):
            if len(toBitString(i, n))>math.floor(math.log2(n)):
                num=i
                break
        if abs(num-n+num) < maxr:
            maxr=abs(num-n+num)
            maxn=n
            maxnum=num
                
    print(maxn)
    ns.append(maxn)
    thrs.append(maxnum)

print(ns)
print(thrs)
"""

"""
maxn=0
maxr=4000
maxnum=0
for n in range(2048, 4096):
    num=0
    for i in range(n):
        if len(toBitString(i, n))>math.floor(math.log2(n)):
            num=i
            break
    if abs(num-n+num) < maxr:
        maxr=abs(num-n+num)
        maxn=n
        maxnum=num
            
            
print(maxn, maxnum)
"""

"""
ns=[21, 43, 85, 171, 341, 683, 1365, 2731, 5461, 10923]
thrs=[11, 21, 43, 85, 171, 341, 683, 1365, 2731, 5461]

best=0
minR=-1

e=0
for i in range(len(ns)):
    res=0
    files = os.listdir("testni-nabor")
    for path in files:
        #print("testni-nabor/" + path)
        try:
            f = open("testni-nabor/" + path, "r")
            vhod = f.read()

            izhod, R = naloga2(vhod, 0, ns[i], thrs[i])
            dek, Rd = naloga2(izhod, 1, ns[i], thrs[i])
            res+=Rd
            print(res)
        except:
            e=1 

    print(res)

    if minR==-1 or res<minR:
        minR=res
        best=ns[i]

print("----")
print(best)
"""

"""
o = open("out", "wb")
buf=""
writeToFile("0000011111010011", o, buf)

o.close()
i = open("out", "rb")

fromBitString(i, 10)
"""
