## get GC content
def checkGCcontent(str seq):
    cdef int total_bases
    cdef double gc_frac, gc_bases
    total_bases = len(seq)
    gc_bases = len([x for x in seq if x == 'C' or x == 'G'])
    gc_frac = gc_bases/total_bases
    return gc_frac

## search matching string and return the position, hamming dist, string
cdef int ham_dist(str s1, str s2) except *:
    if len(s1) != len(s2):
        raise ValueError("Undefined")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def search_min_dist(str source, str search):

    cdef int l, min_dist, index, d
    l = len(search)
    index = 0
    min_dist = l
    cdef str min_substring = source[:l]

    for i in range(len(source)-l+1):
        d = ham_dist(search, source[i:i+l])
        if d<min_dist:
            min_dist = d
            index = i
            min_substring = source[i:i+l]
    return min_dist
