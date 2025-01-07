import numpy as np
a=np.array([1,2,3,4],dtype="int64")
b=np.array([2,3,4,5],dtype="int8")
print(a.dtype,b.dtype)
c=a+b
print(c,c.dtype)
