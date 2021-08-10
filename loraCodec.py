# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 16:41:47 2021

@author: du fanping
"""

#import scipy as sp
import numpy as np
from codecFunctions import whiteningLora
from codecFunctions import crcLora
from codecFunctions import hammingLora
from codecFunctions import xLora
from codecFunctions import grayCode

#system parameters
headH = [
    [1,1,1,1,0,0,0,0,0,0,0,0],
    [1,0,0,0,1,1,1,0,0,0,0,1],
    [0,1,0,0,1,0,0,1,1,0,1,0],
    [0,0,1,0,0,1,0,1,0,1,1,1],
    [0,0,0,1,0,0,1,0,1,1,1,1]
    ]

#'0xB8'
whiteningPoly = [1,0,1,1,1,0,0,0]

#'0x1021'
crcPoly = [0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1]

hammingH = [[1,1,0,1],[1,0,1,1],[1,1,1,0],[0,1,1,1]]

#configure paramters
SF = 12

reduceMode = 'off'
headMode = 'explicit'

nPayload = 255
CR = 1
crcEnb = 1

nCodewordInHead = 0
nSymbolInHead = 0
nHeadCodeword = 0

if SF<7:
    headMode = 'implicit'

if headMode=='explicit':
    nCodewordInHead = SF - 2
    nSymbolInHead = 8
    nHeadCodeword = 5

if SF<11:
    reduceMode = 'off'

if reduceMode=='on':
    SF = SF - 2
    
nCodeword = nPayload*2 + crcEnb*4 - (nCodewordInHead - nHeadCodeword)
nBlock = int(np.ceil(nCodeword/SF))
nSymbol = nBlock*(4+CR) + nSymbolInHead

#payload for test
txBytes = np.random.randint(0,256,nPayload)
txBits = np.zeros((nPayload,8))
for i in range(0,nPayload):
    txBits[i][:] = list(map(int,bin(txBytes[i])[2:].zfill(8)))

#encode
#whitening    
whiteningData = whiteningLora(txBits,whiteningPoly)
whiteningNibble = np.zeros((nPayload*2,4))
for i in range(0,nPayload):
    whiteningNibble[2*i][:] = whiteningData[i][4:8]
    whiteningNibble[2*i+1][:] = whiteningData[i][0:4]
hammingNibble = whiteningNibble
    
#head encode
if headMode=='explicit':
    headBit = list(map(int,bin(nPayload*16+CR*2+crcEnb)[2:].zfill(12)))
    headCheck = list(np.dot(headH,headBit)%2)
    headNibble = np.array([headBit[0:4],headBit[4:8],headBit[8:12],[0,0,0,headCheck[0]],headCheck[1:5]])
hammingNibble = np.vstack((headNibble,hammingNibble))

#crc    
if crcEnb:
    crc = crcLora(txBits,crcPoly)
    crcNibble = np.array([crc[12:16],crc[8:12],crc[4:8],crc[0:4]])
hammingNibble = np.vstack((hammingNibble,crcNibble))
    
#hamming encode
hammingCodeword = hammingLora(hammingNibble,hammingH,CR,'enc')
    
#filling
if nCodeword%SF!=0:
    hammingCodeword = np.vstack((hammingCodeword,np.zeros((SF-nCodeword%SF,8))))

#interleave & symbol gray2bin
if reduceMode=='on':
    txSymbol = np.zeros((nSymbol,SF+2))
else:
    txSymbol = np.zeros((nSymbol,SF))

if headMode=='explicit':
    symbol = xLora(hammingCodeword[0:nCodewordInHead][:],nCodewordInHead,4,1,'x')
    for i in range(0,nSymbolInHead):
        txSymbol[i][0:nCodewordInHead] = grayCode(symbol[i][:],'g2b')
    
for i in range(0,nBlock):
    symbol = xLora(hammingCodeword[nCodewordInHead+i*SF:nCodewordInHead+(i+1)*SF][:],SF,CR,0,'x')
    for j in range(0,4+CR):
        txSymbol[nSymbolInHead+i*(4+CR)+j][0:SF] = grayCode(symbol[j][:],'g2b')

#channel
rxSymbol = txSymbol

#clear encode variable
del txSymbol
del symbol
del hammingCodeword
del hammingNibble
if crcEnb:
    del crcNibble
    del crc
    del txBits
if headMode=='explicit':
    del headNibble
    del headCheck
    del headBit
del whiteningNibble
del whiteningData
del txBytes

#decode
hammingNibble = np.zeros((nCodewordInHead + nBlock*SF,4))
#head
if headMode=='explicit':
    symbol = np.zeros((nSymbolInHead,nCodewordInHead))
    for i in range(0,nSymbolInHead):
        symbol[i][:] = grayCode(rxSymbol[i][0:nCodewordInHead],'b2g')
    hammingCodeword = xLora(symbol,nCodewordInHead,4,1,'o')
    for i in range(0,nCodewordInHead):
        hammingNibble[i][:] = hammingLora(hammingCodeword[i][:],hammingH,4,'dec')
    headNibble = hammingNibble[0:nHeadCodeword][:]
    headBit = np.array(headNibble[0:3][:]).flatten()
    headCheck = np.dot(headH,headBit)%2
    if np.sum(abs(headCheck-np.hstack((headNibble[3][3],headNibble[4][:]))))!=0:
        print('Head Error!')

#payload
#symbol bin2gray 
symbol = np.zeros((nSymbol-nSymbolInHead,SF))
for i in range(0,nSymbol-nSymbolInHead):
    symbol[i][:] = grayCode(rxSymbol[i+nSymbolInHead][0:SF],'b2g')

#deinterleave
for i in range(0,nBlock):
    hammingCodeword = xLora(symbol[i*(4+CR):(i+1)*(4+CR)][:],SF,CR,0,'o')
    #hamming decode
    for j in range(0,SF):
        hammingNibble[nCodewordInHead+i*SF+j][:] = hammingLora(hammingCodeword[j][:],hammingH,CR,'dec')

#whitening
whiteningData = np.zeros((nPayload,8))        
for i in range(0,nPayload):
    whiteningData[i][:] = np.hstack((hammingNibble[nHeadCodeword+2*i+1][:],hammingNibble[nHeadCodeword+2*i][:]))    
rxBits = whiteningLora(whiteningData,whiteningPoly)

#crc
if crcEnb:
    crc = crcLora(rxBits,crcPoly)
    crcNibble = hammingNibble[nHeadCodeword+nPayload*2:nHeadCodeword+nPayload*2+4][:]
    rxCrc = np.hstack((crcNibble[3][:],crcNibble[2][:],crcNibble[1][:],crcNibble[0][:]))
    if np.sum(abs(crc-rxCrc))!=0:
        print('Crc Error!')
elif np.sum(abs(rxBits-txBits))!=0:
    print('Codec Error!')
    
