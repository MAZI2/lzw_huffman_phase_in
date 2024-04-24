import math
import os
import subprocess
import sys
import numpy as np
import heapq

dictLen = 4096
blockSize = 100

# increase recursion limit because of the Huffman trees
sys.setrecursionlimit(10000)

# Huffman tree - code lengths
class node:
    def __init__(self, freq, symbol, left=None, right=None):
        self.freq = freq
        self.symbol = symbol
        self.left = left
        self.right = right
        self.huff = ''

    def __lt__(self, nxt):
        return self.freq < nxt.freq

def printNodes(node, lengths, val=''):
    newVal = val + str(node.huff)
    if(node.left):
        printNodes(node.left, lengths, newVal)
    if(node.right):
        printNodes(node.right, lengths, newVal)
    if(not node.left and not node.right):
        lengths[node.symbol]=len(newVal)

def hufLengths(freq):
    nodes = []

    for x in range(len(freq)):
        heapq.heappush(nodes, node(freq[x], x))

    while len(nodes) > 1:
        left = heapq.heappop(nodes)
        right = heapq.heappop(nodes)

        left.huff = 0
        right.huff = 1

        newNode = node(left.freq+right.freq, left.symbol+right.symbol, left, right)

        heapq.heappush(nodes, newNode)

    lengths=[0 for i in range(len(freq))]
    printNodes(nodes[0], lengths)
    return lengths

# construct canonical Huffman code from given lengths
def huffman(length):
    syms=len(length)
    count=[]
    for l in range(16):
        count.append(0)
    for sym in range(syms):
        count[length[sym]]+=1
    if count[0] == syms:
        return 0

    count[0] = 0
    code = 0
    next_code = []
    next_code.append(0)
    for bits in range(1, max(length)+1):
        code = (code + count[bits-1]) << 1
        next_code.append(code)

    codes = dict()
    for n in range(syms):
        l = length[n]
        if l != 0:
            codes[n]=(next_code[l])
            next_code[l]+=1


    for i in range(syms):
        if length[i] != 0:
            form = '{0:0' + str(length[i]) + 'b}'
            codes[i]=form.format(codes[i])
    return codes, count

# decode from canonical Huffman code 
# returns index of dictionary element
def decode(s, lenFreqs, pos=0):
    code = 0
    first = 0
    index = 0
    p = pos

    for l in range(1, 16):
        num = lenFreqs[l]
        code = code << 1
        code += int(s[p])
        p+=1
        first = first << 1
        first += num
        index += num
        if code < first:
            return index+(code-first), p

# decode Huffman code for code lengths of Huffman code for the encoded block data
# decode the lengths for encoded Huffman code for data
# decode and return block data
def decodeBlock(HLIT, huffman2, buf, file, max_acc):
    order = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]

    codeCodeLengths = [0 for i in range(19)]
    i=0
    for o in order:
        codeCodeLengths[o]=huffman2[i]
        i+=1

    # count ... frequencies of lengths
    cclHuff, count = huffman(codeCodeLengths)

    symbols = sorted(cclHuff, key=lambda key: len(cclHuff[key]))
    lengths_final = []

    pos = 0
    prev = -1
    first = -1

    debugBuf = buf
    # control the amount of bytes read
    c = 0

    # decode lengths for the data Huffman code (Huffman + RLE-ish decode)
    for i in range(HLIT):
        if c>=HLIT:
            break
        if len(buf) < 8:
            s=''.join(format(byte, '08b') for byte in file.read(1))
            debugBuf+=s
            buf+=s

        si, p = decode(buf, count, pos)
        buf = buf[p:]

        if len(buf) < 8:
            s=''.join(format(byte, '08b') for byte in file.read(1))
            debugBuf+=s
            buf+=s

        sym = symbols[si]
        if first == -1:
            first = sym

        extra = 0
        if sym == 16:
            c+=1
            extra = int(buf[:2], 2)+3
            buf = buf[2:]
            for j in range(extra):
                lengths_final.append(prev)
        elif sym == 17:
            c+=1
            extra = int(buf[:3], 2)+3
            buf = buf[3:]
            for j in range(extra):
                lengths_final.append(first)
        elif sym == 18:
            c+=1
            extra = int(buf[:7], 2)+11
            buf = buf[7:]
            for j in range(extra):
                lengths_final.append(first)
        else:
            lengths_final.append(sym)

        c+=1
        prev = sym

#    print(debugBuf)

#    print(lengths_final)

    # create data Huffman code from decoded lengths and decode block data
    huffData, count = huffman(lengths_final)

#    print(huffData)
    data_symbols = sorted(huffData, key=lambda key: len(huffData[key]))
#    print(data_symbols)
    block_decoded = []

    pos = 0
    debugbuf2=buf
#    print("TERM", huffData[max_acc+1])

    # decode data until the terminator symbol is decoded
    # terminator sybol = max access index + 1 (appended to block data)
    while True:
        if len(buf) < max(lengths_final):
            s=''.join(format(byte, '08b') for byte in file.read(1))
            buf+=s
            debugbuf2+=s
        si, p = decode(buf, count, pos)
        buf=buf[p:]
        sym = data_symbols[si]
        if sym == max_acc+1:
#            print("End of block reached")
            break
        block_decoded.append(sym)

#   print(debugbuf2)
#   print(buf)

    return block_decoded, buf

# encode dictinary accesses (cl)
# RLE-ish encode lengths for data Huffman code
# encode the encoded lengths (ccl)
# return HLIT = number of ccl, ccl Huffman code in standard format, ccl bit string, encoded block data bit string
def encodeBlock(dict_acc, freqs):
    # get lengths of Huffman code from dict_acc frequencies (freqs)
    filtered_freqs = [x for x in freqs if x != 0]
    norm_lengths = hufLengths(filtered_freqs)

    lengths = freqs
    for i in range(len(freqs)):
        if freqs[i] > 0:
            lengths[i]=norm_lengths.pop(0)

#    print(lengths)

    block_data = dict_acc

    # RLE-ish encoding of lengths
    freqs2=[0 for i in range(19)]
    first = lengths[0]
    freqs2[first]+=1
    seq = []
    seq.append([first, 0])
    repNum = 0
    firNum = 0
    prev = first
    for l in lengths[1:]:
        if l == prev and l != first and repNum<=5:
            repNum+=1
            prev = l
            continue
        elif repNum>=3:
            seq.append([16, repNum])
            freqs2[16]+=1
            repNum=0
        else:
            for i in range(repNum):
                seq.append([prev, 0])
                freqs2[prev]+=1
            repNum=0

        if l == first and firNum<138:
            firNum+=1
            prev = l
            continue
        elif firNum>=11:
            seq.append([18, firNum])
            freqs2[18]+=1
            firNum=0
        elif firNum>=3:
            seq.append([17, firNum])
            freqs2[17]+=1
            firNum=0
        else:
            for i in range(firNum):
                seq.append([prev, 0])
                freqs2[prev]+=1
            firNum=0

        seq.append([l, 0])
        freqs2[l]+=1

        prev = l

    if repNum>=3:
        seq.append([16, repNum])
        freqs2[16]+=1
    else:
        for i in range(repNum):
            seq.append([prev, 0])
            freqs2[prev]+=1

    if firNum>=11:
        seq.append([18, firNum])
        freqs2[18]+=1
    elif firNum>=3:
        seq.append([17, firNum])
        freqs2[17]+=1
    else:
        for i in range(firNum):
            seq.append([prev, 0])
            freqs2[prev]+=1

#    print(seq)

    # Huffman code for code lengths (ccl)
    lengths2 = hufLengths(freqs2)

    for l in range(len(lengths2)):
        if freqs2[l]==0:
            lengths2[l]=0

#    print(lengths2)

    huff, _ = huffman(lengths2)

#    print(huff)

    # encode the RLE-ish sequence sequence
    encodedLengths = []
    for s in seq:
        encodedLengths.append(huff[s[0]])
        if s[0]==16:
            encodedLengths.append('{0:02b}'.format(s[1]-3))
        elif s[0]==17:
            encodedLengths.append('{0:03b}'.format(s[1]-3))
        elif s[0]==18:
            encodedLengths.append('{0:07b}'.format(s[1]-11))

    # format in standard order
    order = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]
    encodedHuffman = []
    for o in order:
        encodedHuffman.append(lengths2[o])

#    print(encodedHuffman)
#    print(encodedLengths)

    # encode block data
    # construct data Huffman tree
    huffPrim, _ = huffman(lengths)

#    print(huffPrim)

    block_encoded = []
    for c in block_data:
        block_encoded.append(huffPrim[c])

    block_encoded="".join(block_encoded)

    # format HLIT
    HLIT = '{0:08b}'.format(len(encodedLengths))

#    print(HLIT)

    # format the ccl Huffman tree as 3bit number sequence
    huffman2 = []
    for i in encodedHuffman:
        huffman2.append('{0:03b}'.format(i))
    huffman2 = "".join(huffman2)

#    print(huffman2)

    encodedLengths = "".join(encodedLengths)

#    print(encodedLengths)
#    print(block_encoded)

    return HLIT, huffman2, encodedLengths, block_encoded

# read HLIT, ccl Huffman tree sequence, max_acc (for terminator symbol)
def extractHeader(file, buf):
    #outputs numbers
    izhod = []
    max_acc = 0
    for i in range(9):
        b=file.read(1)
        temp = buf
        buf+=''.join(format(byte, '08b') for byte in b)
        if buf == temp:
            return "", "", "", 0, 0

    HLIT = int(buf[:8], 2)
    buf = buf[8:]

    huffman2 = []
    for i in range(19):
        huffman2.append(int(buf[:3], 2))
        buf = buf[3:]

    while len(buf)<12:
        buf+=''.join(format(byte, '08b') for byte in file.read(1))

    max_acc = int(buf[:12], 2)
    buf = buf[12:]

    return HLIT, huffman2, buf, max_acc, 1

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

# assemble the whole block
def constructBlock(dict_acc, freqs, bout, max_acc, buf):
    dict_acc.append(max_acc+1)
    freqs.append(1)
    HLIT, huffman2, encodedLengths, block_encoded = encodeBlock(dict_acc, freqs)
    term = format(max_acc, '012b')

    buf = writeToFile(HLIT, bout, buf)
    buf = writeToFile(huffman2, bout, buf)
    buf = writeToFile(term, bout, buf)
    buf = writeToFile(encodedLengths, bout, buf)
    buf = writeToFile(block_encoded, bout, buf)
    return buf

# dissasemle the encoded block
def deconstructBlock(blin, fout, buf):
    HLITd, huffman2d, buf, max_acc, hasNext = extractHeader(blin, buf)
    if hasNext == 0:
        return

    decoded_block, buf = decodeBlock(HLITd, huffman2d, buf, blin, max_acc)
    # LZW decode from decoded dictionary access sequence
    dekodiranje(decoded_block, dictLen, "d", fout)

    # recursive call to decode all blocks
    deconstructBlock(blin, fout, buf)
    return

# -------------
# LZW ... TODO: translate to English

# prefill the LZW dictionary
def loadDict(typ, f, dec):
    preDict=[]
    header=[]
    
    slovar=dict()

    # encoding / decoding
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

# LZW decompression, writes output to file
def dekodiranje(vhod: list, n, typ, out) -> list:
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

    for c in izhod:
        out.write(ord(c).to_bytes(1, byteorder='big'))

# LZW compression, repeats recursively for blocks
def kodiranje(vhod: list, n, typ, file, buf) -> list:
    slovar, header = loadDict(typ, vhod, 0)
    bcounter = 0 

    index = len(slovar)
    #if len(header) != 0:
        
    accs = []
    max_acc = -1;
    N = '' 
    while b:=vhod.read(1):
        z = chr(int.from_bytes(b, byteorder='big'))
        Nz = N + z 
        if Nz in slovar:
            N = Nz
        else:
            accs.append(slovar[N])
            if(slovar[N] > max_acc):
                max_acc = slovar[N]

            if len(slovar) + 1 < dictLen:
                slovar[Nz] = index 
                index += 1;
            N = z

        bcounter+=1
        if bcounter==blockSize:
            break

    accs.append(slovar[N])
    if(slovar[N] > max_acc):
        max_acc = slovar[N]

    freqs = [0 for i in range(max_acc+1)]
    for a in accs:
        freqs[a]+=1

    buf = constructBlock(accs, freqs, file, max_acc, buf)

    # empty the buffer (extend to 8 bits and write) if reached end of file
    if bcounter<blockSize:
        writeToFile("", file, buf)
        return
    kodiranje(vhod, n, typ, file, buf)

    return accs, freqs, max_acc

# driver function
def naloga2_tekma(dat_vhod: str, dat_izhod: str, nacin: int) -> float:

    # unneeded (file type if you want different dictionary preffils)
    typ = "d"

    # encode
    buf = ""
    inp = open(dat_vhod, "rb")
    bout = open("code.lzw", "wb")

    accs, freqs, max_acc = kodiranje(inp, dictLen, typ, bout, buf)
    bout.close()

    # decode
    buf = ""
    fin = open("code.lzw", "rb")
    fout = open("decoded_file", "wb")

    deconstructBlock(fin, fout, buf)
    fout.close()

    return 0

files = os.listdir("testfiles")

for path in files:
    abs_path = "testfiles/" + path

    print(path)
    naloga2_tekma(abs_path, os.path.splitext(path)[0] + ".lzw", 0)
