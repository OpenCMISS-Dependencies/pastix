from Projects.PaStiX import PaStiX
from Projects.Hips import Hips
def printCompilationNames(p):
    """ print compilation keywords"""
    
    out = ''
    for l in p.listCompilationNames():
        out += l.pop(0) + u'\t\t => '
        if l:
            w = l.pop(0)
        else:
            w = None
        while not w == None:
            out += w
            if l:
                w = l.pop(0)
                out += u' : '
            else:
                w = None
        out += '\n'
    print out


print """
=== PaStiX Compilation options ===
"""

p = PaStiX()
printCompilationNames(p)
print """
=== Hips Compilation options ===
"""
p = Hips()
printCompilationNames(p)

