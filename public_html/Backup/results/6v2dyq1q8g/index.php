Traceback (most recent call last):
  File "./script.py", line 443, in <module>
    start_sampling(seq)
  File "./script.py", line 368, in start_sampling
    smp = sampling(sequence, seq.RNA_entry, seq.fixedCG_entry)
  File "./script.py", line 87, in __init__
    self.wild_entry = getbpp(fname)
  File "./script.py", line 139, in getbpp
    f = open(fname+'_dp.ps','r')
IOError: [Errno 2] No such file or directory: '92-52_dp.ps'