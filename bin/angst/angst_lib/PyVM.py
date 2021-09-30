import pdb
import os

_proc_status = '/proc/%d/status' % os.getpid()
_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}
_total_memory = 0.0
_last_memory = 0.0

def _VmB(VmKey):
    '''private method'''
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to MB
    return (float(v[1]) * _scale[v[2]]) / (1024.0*1024.0)


def MemoryUpdate(print_msg=None,str_return=False):
    '''print memory usage stats in MB.
    '''
    global _total_memory, _last_memory

    _last_memory = _total_memory
    _total_memory = _VmB('VmSize:')

    if print_msg is not None:
        mem_diff = _total_memory - _last_memory
        m1 = mem_diff
        m2 = _total_memory
        print_str = '%3.3f\t%3.3f\t' % (m2, m1)
        print_str += print_msg
        if str_return:
            return print_str + "\n"
        else:
            print print_str



