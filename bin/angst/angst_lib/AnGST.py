#!/usr/bin/python
# AnGST

# python libraries
import sys
import time
import pdb
import PyVM

from AnGSTHelper import RunAnGST
from AnGSTInput import input_obj

# initiate variable for measuring running time and memory consumtion
start_time = time.time()
mem_str = []
mem_str.append(PyVM.MemoryUpdate("init",'return'))

# load inputs #
print "* read input"
input_info = input_obj(sys.argv[1])
input_dict = {}
input_dict['input_info'] = input_info
input_dict['mem_str'] = mem_str
input_dict['write_out'] = False
angst_inputs = [input_dict]

# run angst one more time at that scaling
input_dict['write_out'] = True
input_dict['start_time'] = start_time
RunAnGST(input_dict)


