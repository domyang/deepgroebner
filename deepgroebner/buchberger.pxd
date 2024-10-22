from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.utility cimport pair

cdef extern from "buchberger.cpp":
    pass

cdef extern from "buchberger.h":
    cdef struct Numpy1DArray:
        int* ptr
        int size

    cdef cppclass LeadMonomialsEnv:
        LeadMonomialsEnv() except +
        LeadMonomialsEnv(vector[Numpy1DArray], int*, bint, bint, int, double*, string, int, double) except +
        LeadMonomialsEnv(const LeadMonomialsEnv&)
        void reset(int*, double*, double)
        pair[double,bint] step(int)
        void seed(int)
        int select(int)
        vector[double] get_cost()
        double objective()
        double value(string, double, int)

        vector[int] state
        vector[int] basis
        vector[int] pairs
        int cols

