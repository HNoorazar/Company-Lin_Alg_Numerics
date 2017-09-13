import numpy as np
def createTent(center):
    assert center >= 0.0 and center <= 1.0
    def tent(x):
        if x >= 0 and x < center:
            return 1.0 / center
        elif x >= center and x <= 1.0:
            return -1.0 / (1.0 - center)
        else:
            return 0
    return tent

tent_func = createTent(.3)
print tent_func(.1)
print tent_func(.2)
"""
v = np.arange(0,6)/10.
print "sin(x)= ", sin(v)
print "type", type(v)
print "tent_func(v)", tent_func(v)
"""

class MyClass(object):
     i = 123
     def __init__(self):
         self.i = 345

a = MyClass()
print a.i
345
print MyClass.i
123


class SomeClass:
    def __init__(self, y):
        self.y = y
 
    def foo(self, x):
        return self.y + x	
a = SomeClass(3)
b = SomeClass(4)
print(a.foo(5))
print(b.foo(5))

class SomeOtherClass:
   def __init__(self, x):
        self.x = x

   @classmethod
   def foo(cls, x):
       return cls(x)

c = SomeOtherClass.foo(5)
print c


a = 1
if a == 0x000f:
   print "yes"
else:
   print "no"
print 0x000f

