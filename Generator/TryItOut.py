#!/usr/bin/env python

from BlasTensorContraction import Tensor

# CLASS 2
T1 = Tensor(["a", "b", "c"])
T2 = Tensor(["b", "c"])
print T1 * T2

T1 = Tensor(["a", "b", "c"])
T2 = Tensor(["a", "b"])
print T1 * T2

T1 = Tensor(["a", "b", "c"])
T2 = Tensor(["a", "c"])
print T1 * T2

# CLASS 3
T1 = Tensor(["a", "b", "c"])
T2 = Tensor(["n", "b", "c"])
print T1 * T2

T1 = Tensor(["a", "b", "c"])
T2 = Tensor(["c", "a", "n"])
print T1 * T2

# CCD and GR
T1 = Tensor(["i", "a", "j", "b"])
T2 = Tensor(["i", "c", "j", "d"])
print T1 * T2

T1 = Tensor(["i", "a", "j", "b"])
T2 = Tensor(["j", "c", "i", "d"])
print T1 * T2

# test
T1 = Tensor(["i", "j", "k", "l", "m", "n"])
T2 = Tensor(["j", "k", "l", "z", "v", "i"])
print T1 * T2

# test 1
T1 = Tensor(["i", "j", "k", "l", "m", "n"])                          
T2 = Tensor(["j", "k", "l", "z", "v", "n"])                                                                                              
print T1 * T2

# test 2
T1 = Tensor(["i", "j", "k", "l", "m", "n"])
T2 = Tensor(["s", "n", "j", "z", "v", "i"])
print T1 * T2

# test 3
T1 = Tensor(["i", "j", "k", "l", "m", "n"])                                                                                              
T2 = Tensor(["k", "n", "j", "z", "l", "i"])                                                                                              
print T1 * T2

# test 4
T1 = Tensor(["i", "j", "k", "l", "m", "n"])
T2 = Tensor(["k", "n", "j", "m", "l", "i"])
print T1 * T2

# test 5
T1 = Tensor(["i", "j", "k", "l", "m", "n"])
T2 = Tensor(["k", "n", "j", "l", "i"])
print T1 * T2
