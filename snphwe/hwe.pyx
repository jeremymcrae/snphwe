# distutils: sources = src/snp_hwe.cpp
# distutils: include_dirs = src/

cdef extern from "snp_hwe.h":
    # double SNPHWE(int, int, int) except +
    double SNPHWE(int, int, int)

def snphwe(obs_hets, obs_hom1, obs_hom2):
    return SNPHWE(obs_hets, obs_hom1, obs_hom2)
