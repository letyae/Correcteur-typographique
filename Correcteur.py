#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np




T = np.zeros((27, 27), float)


code = 0
with open("proust.txt", "r") as f:
        line = f.readlines()
        for s in line:
            for c in s:
                if c.isalpha():
                    T[code, ord(c)-96] += 1
                    code = ord(c)-96
                else:
                    T[code, 0] += 1
                    code = 0
                        
for i in range(T.shape[0]):
    T[i] /= T[i].sum()
    

    
    


def encode(c):
    if c.isalpha():
        return(ord(c)-96)
    else:
        return(0)
    
def encode_s(s):
    l = []
    for c in s:
        l.append(encode(c))
    return(l)

def decode_s(l):
    s = ""
    for e in l:
        if e == 0:
            s += " "
        else:
            s += chr(e+96)
    return(s)



def p_stat(T):
    p = np.random.randn(len(T))
    p = p@np.linalg.matrix_power(T,100)
    return(p/p.sum())

def p_markov(s,T):
    p_s = p_stat(T)
    l = encode_s(s)
    p = p_s[l[0]]
    for i in range(len(l)-1):
        p = p * T[l[i], l[i+1]]
    return(p)




        
s = "au revoir ma chere"
print(p_markov(s, T))
s = "au revoir ma chrte"
print(p_markov(s, T))
s = "ay revoir mz chere"
print(p_markov(s, T)) 




d = {}

d[encode(' ')] = encode_s(' ')
d[encode('a')] = encode_s('z')
d[encode('b')] = encode_s('vn')
d[encode('c')] = encode_s('xv')
d[encode('d')] = encode_s('sf')
d[encode('e')] = encode_s('zr')
d[encode('f')] = encode_s('dg')
d[encode('g')] = encode_s('fh')
d[encode('h')] = encode_s('gj')
d[encode('i')] = encode_s('uo')
d[encode('j')] = encode_s('hk')
d[encode('k')] = encode_s('jl')
d[encode('l')] = encode_s('km')
d[encode('m')] = encode_s('l')
d[encode('n')] = encode_s('b')
d[encode('o')] = encode_s('ip')
d[encode('p')] = encode_s('o')
d[encode('q')] = encode_s('s')
d[encode('r')] = encode_s('et')
d[encode('s')] = encode_s('qd')
d[encode('t')] = encode_s('ry')
d[encode('u')] = encode_s('yi')
d[encode('v')] = encode_s('cb')
d[encode('w')] = encode_s('x')
d[encode('x')] = encode_s('wc')
d[encode('y')] = encode_s('tu')
d[encode('z')] = encode_s('ae')




p_err = 0.1

LT = []
for i in range(len(T)):
    LT.append(np.zeros(T.shape))
    LT[-1][:,i] += T[:,i]*(1-p_err)
    for c in d[i]:
        LT[-1][:,c] += T[:,c]*(p_err)/len(d[c])

        

LT = np.array(LT)
    




def p_hmm(s,LT):
    T = LT.sum(axis = 0)
    p_i = p_stat(T)
    p_t = np.ones(len(T))
    l = encode_s(s)
    for i in range(len(l)):
        p_i = p_i @ LT[l[i]]
    return(p_i@p_t)




s = "au revoir ma chere"
print(p_hmm(s, LT))
s = "au revoir ma chrte"
print(p_hmm(s, LT))
s = "ay revoir mz chere"
print(p_hmm(s, LT)) 
    




def viterbi(s, LT):
    p = p_stat(LT.sum(axis=0))

    p_max = [p]
    c_max = [[[]]*27]


    observations = encode_s(s)

    for o in observations[:]:
        p_new = np.array(p, float)
        l = []
    
        for i in range(LT.shape[0]):
            p_0 = p_max[-1]*LT[o,:,i]
            p_new[i] = np.max(p_0)
            l.append(np.argmax(p_0))
    
        p_max.append(p_new)
    
        c_maxold = list(c_max[-1])
        c_maxnew = list(c_max[-1])
    
        for i in range(LT.shape[0]):
            c_maxnew[i] = c_maxold[l[i]] + [i]
    
        c_max.append(c_maxnew)
        
    return(p_max, c_max)



s = "dabs la pezieoz lrs vhiens cputebt zpres les miytonq ey les bkrs sont tebsees"


p_max,c_max= viterbi(s, LT)

l = c_max[-1][np.argmax(p_max[-1])]
s2 = decode_s(l)
print(s2)

