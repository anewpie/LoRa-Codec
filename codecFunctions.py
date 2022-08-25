# -*- coding: utf-8 -*-
"""
Created on Aug 11 2021
@author: anewpie
"""

import numpy as np

#LFSR:'B8'
def whiteningLora(inBits,whiteningPoly):
    
    regs = np.ones(8)
    outBits = np.zeros((len(inBits),8))
    
    for i in range(0,len(inBits)):
        outBits[i][:] = (inBits[i][:]+regs)%2
        regs = np.hstack((regs[1:8],np.dot(regs,whiteningPoly)%2))
        
    return outBits

#CCITT16:'1021'
def crcLora(inBits,genPoly):

    crc = np.zeros(16)

    for i in range(0,len(inBits)):
        for j in range(0,8):
            crc = (np.hstack((crc[1:16],inBits[i][j]))+crc[0]*np.array(genPoly))%2
        
    return crc

#hamming(8,4)
def hammingLora(inBits,H,CR,mode):
    
    if mode=='enc':
        outBits = np.zeros((len(inBits),8))
        for i in range(0,len(inBits)):
            outBits[i][:] = np.dot(np.vstack((np.array(H),np.eye(4))),inBits[i][:])%2
    elif mode=='dec':
        outBits = inBits[CR:4+CR]
        syndrome = np.dot(np.hstack((np.eye(CR),np.array(H[4-CR:4][0:4]))),inBits[:])%2
        if CR>2:
            check = ((np.outer(syndrome,np.ones(4))+np.array(H[4-CR:4][0:4]))%2).T
            for i in range(0,4):
                if np.sum(check[i][:])==0:
                    outBits[i] = (outBits[i]+1)%2
                
    return outBits

#block x
def xLora(inBits,SF,CR,mode):
    
    if mode=='x':
        outBits = np.flipud(np.fliplr(inBits.T))
        #special for CR=1
        if CR==1:
            outBits[4][:] = (outBits[4][:] + outBits[3][:])%2
        for i in range(0,4+CR):
            outBits[i][:] = np.hstack((outBits[i][SF-i%SF:SF],outBits[i][0:SF-i%SF]))
    elif mode=='o':
        for i in range(0,4+CR):
            inBits[i][:] = np.hstack((inBits[i][i%SF:SF],inBits[i][0:i%SF]))
        #special for CR=1
        if CR==1:
            inBits[4][:] = (inBits[4][:] + inBits[3][:])%2
        outBits = np.flipud(np.fliplr(inBits.T))
    
    return outBits

#gray
def grayCode(inBits,mode):
        
    if mode=='b2g':
        nIter = 1
    elif mode=='g2b':
        nIter = int(np.ceil(np.log2(len(inBits))))

    for i in range(0,nIter):
        inBits = (inBits+np.hstack((np.zeros(2**i),inBits[0:len(inBits)-2**i])))%2

    return inBits    

