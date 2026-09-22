"""Numerical check of Lemma 1 (staggered trapezoids) of T1 v0.4: (a) K=3 exact formula, (b) sufficient condition for K<=5, (c) K=4 upper bound.
Run: python3 sims/check_lemma1_trapezoid.py   (needs check_thm1_support.py in the same folder)"""
import numpy as np, itertools, sys
from math import gcd
from functools import reduce
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # check_thm1_support.py lives beside this file
from check_thm1_support import identified_set
def nullity(entries, T): return identified_set(entries, T)[0].shape[0]
def D(entries): return reduce(gcd, [abs(a-b) for a in entries for b in entries if a!=b])

print("(a) K=3: exact formula nullity = max(d, (e2-e1) - (T-e3)) ?")
bad=0; n=0
for entries in itertools.combinations(range(1,16), 3):
    entries=list(entries)
    if entries[0]!=1: continue
    d=D(entries); D2=entries[1]-entries[0]
    for T in range(entries[-1]+1, entries[-1]+18):
        w=T-entries[-1]; pred=max(d, D2-w); got=nullity(entries,T); n+=1
        if pred!=got:
            bad+=1
            if bad<=8: print("   mismatch:", entries, "T=",T, "w=",w, "d=",d, "pred",pred,"got",got)
print(f"   checked {n}, mismatches {bad}")

print("\n(b) all K in {2,3,4,5}: sufficient condition  T - e_K >= (e2-e1) - d  =>  nullity = d ?")
bad=0; n=0; bad_rows=[]
for K in [2,3,4,5]:
    for entries in itertools.combinations(range(1,14), K):
        entries=list(entries)
        if entries[0]!=1: continue
        d=D(entries); D2=entries[1]-entries[0]
        for T in range(entries[-1]+1, entries[-1]+14):
            w=T-entries[-1]
            if w >= D2-d:
                n+=1; got=nullity(entries,T)
                if got!=d: bad+=1; bad_rows.append((entries,T,w,d,got))
print(f"   checked {n} (condition holds), violations {bad}"); [print("   ",r) for r in bad_rows[:8]]

print("\n(c) K=4: does  max(d, (e2-e1) - (T-e_K))  still give the nullity? (upper bound? exact?)")
over=under=eq=0; ex=[]
for entries in itertools.combinations(range(1,14), 4):
    entries=list(entries)
    if entries[0]!=1: continue
    d=D(entries); D2=entries[1]-entries[0]
    for T in range(entries[-1]+1, entries[-1]+12):
        w=T-entries[-1]; pred=max(d, D2-w); got=nullity(entries,T)
        if got==pred: eq+=1
        elif got<pred: under+=1; ex.append(('got<pred',entries,T,w,d,pred,got))
        else: over+=1; ex.append(('got>pred',entries,T,w,d,pred,got))
print(f"   equal {eq}, nullity < formula {under}, nullity > formula {over}")
for r in ex[:6]: print("   ", r)
