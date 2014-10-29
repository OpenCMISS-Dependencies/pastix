import re
import sys
in_enum = False
enum_re = re.compile(u".*enum\s*(\w+)\s*{")
enum_field_re = re.compile(u"\s*(\w+)\s*=\s*(\d+)")
enum_end_re = re.compile(u".*}");
enum_remove_re = re.compile("(API_[^_]*)")

enum_list = {}
enum_list2 = {}
removes = {
    u"IPARM_ACCESS"  : u"IPARM_",
    u"DPARM_ACCESS"  : u"DPARM_",
    u"ERR_NUMBERS"   : u"_ERR",
    u"MODULES"       : u"MOD_",
    u"API_TRACEFMT"  : u"API_TRACE_",
    u"API_FLOAT"     : u"API_",
    u"API_BOOLEAN"   : u"API_",
    u"API_ERASE_CSC" : u"API_CSC_",
    u"API_RHS"       : u"API_"
}

replaces = {
    u"IPARM_ACCESS"    : u"IPARM",
    u"DPARM_ACCESS"    : u"DPARM",
    u"ERR_NUMBERS"     : u"ERR",
    u"API_SYM"         : u"SYM",
    u"API_TASK"        : u"TASK",
    u"API_IO"          : u"IO",
    u"API_BOOLEAN"     : u"BOOLEAN",
    u"API_ORDER"       : u"ORDER",
    u"API_TASK_OLD"    : u"TASK_OLD",
    u"API_RHS"         : u"RHS",
    u"API_VERBOSE"     : u"VERBOSE",
    u"API_BIND_MODE"   : u"BIND_MODE",
    u"API_RAF"         : u"RAF",
    u"API_FLOAT"       : u"FLOAT",
    u"API_TRACEFMT"    : u"TRACEFMT",
    u"API_THREAD_MODE" : u"THREAD_MODE",
    u"API_ERASE_CSC"   : u"CSC",
    u"API_FACT"        : u"FACT"
}
for line in open(sys.argv[1]):
    if not in_enum:
        m = enum_re.search(line)
        if (m):
            current_enum = m.group(1)
            enum_list[current_enum] = {}
            enum_list2[current_enum] = {}
            in_enum=True
            if (current_enum in removes.keys()):
                remove = removes[current_enum]
            else:
                m = enum_remove_re.search(current_enum)
                if (m):
                    remove = m.group(1)+"_"
                else:
                    remove = u"API_"
    else:
        m = enum_end_re.search(line)
        if (m):
            in_enum=False
        else:
            m = enum_field_re.search(line)
            if (m):
                name = re.sub(remove, "", m.group(1))
                enum_list[current_enum][name] = m.group(2)
                enum_list2[current_enum][m.group(1)] = m.group(2)

print u'include "pastix_types.pxi"'

#for enum in enum_list.keys():
#    print u"class ENUM_%s:" % enum
#    for enum_field in enum_list[enum].keys():
#        value = enum_list[enum][enum_field]
#        print u"     %-20s = %s" % (enum_field, value)
#    print u""

print u"class ENUM_API:"
for enum in enum_list2.keys():
    for enum_field in enum_list2[enum].keys():
        value = enum_list2[enum][enum_field]
        print u"     %-20s = %s" % (enum_field, value)
    print u""

print u"from libc.stdlib cimport malloc, free"
print u"""
cimport mpi4py.MPI as MPI
from mpi4py.mpi_c cimport *

# PaStiX Exception
class PastixException(Exception):
    pass

class PastixArithmeticException(PastixException):
    pass

"""
print u"cdef class ZPASTIX:"

#for enum in enum_list.keys():
#for enum in enum_list.keys():
#    if enum in replaces.keys():
#        replacement = replaces[enum]
#    else:
#        replacement = enum
#
#    #print u"    %-20s = ENUM_%s" % (replacement, enum)

print u"    %-20s = ENUM_%s" % (u'API', u'API')
print u""
for line in open(sys.argv[2]):
    print line
print u""

print u"cdef class CPASTIX:"
print u"    API = ENUM_API"
print u""
for line in open(sys.argv[2]):
    print line.replace(u"z_",u"c_").replace(u"double complex", "complex").replace("complex128", "complex64")
print u""

print u"cdef class DPASTIX:"
print u"    API = ENUM_API"
print u""
for line in open(sys.argv[2]):
    print line.replace(u"z_",u"d_").replace(u"double complex", "double").replace("complex128", "float64")
print u""

print u"cdef class SPASTIX:"
print u"    API = ENUM_API"
print u""
for line in open(sys.argv[2]):
    print line.replace(u"z_",u"s_").replace(u"double complex", "float").replace("complex128", "float32")
print u""

print u"cdef class PASTIX:"
print u"    API = ENUM_API"
print u""
for line in open(sys.argv[2]):
    print line.replace(u"z_",u"").replace(u"double complex", "pastix_float_t").replace("np.complex128_t", "np_pastix_float_t")
print u""

#for enum in enum_list.keys():
#    print u"del ENUM_%s" % enum
#print u"del ENUM_API"
