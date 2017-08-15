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

v = np.arange(0,6)/10.
print "sin(x)= ", sin(v)
print "type", type(v)
print "tent_func(v)", tent_func(v)