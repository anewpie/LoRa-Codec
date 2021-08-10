# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 16:45:44 2021

@author: du fanping
"""

import numpy as np
from codecFunctions import hammingLora

hammingH = [[1,1,0,1],[1,0,1,1],[1,1,1,0],[0,1,1,1]]

nPayload = 255*2
CR = 3

txBytes = np.random.randint(0,16,nPayload)

txBits = np.zeros((nPayload,4))
for i in range(0,nPayload):
    txBits[i][:] = list(map(int,bin(txBytes[i])[2:].zfill(4)))
    
encode = hammingLora(txBits,hammingH,CR,'enc')

errorBit = np.random.randint(0,4+CR,nPayload)
errorPattern = np.zeros((nPayload,4+CR))
for i in range(0,nPayload): 
    errorPattern[i][:] = list(map(int,bin(2**errorBit[i])[2:].zfill(4+CR)))

rxBits = np.zeros((nPayload,4+CR))
decode = np.zeros((nPayload,4))
for i in range(0,nPayload):
    rxBits[i][:] = (encode[i][4-CR:8] + errorPattern[i][:])%2
    
for i in range(0,nPayload):    
    decode[i][:] = hammingLora(rxBits[i][:],hammingH,CR,'dec')
    
if np.sum(abs(decode-txBits))!=0:
    print('Decode Error!')
    
    