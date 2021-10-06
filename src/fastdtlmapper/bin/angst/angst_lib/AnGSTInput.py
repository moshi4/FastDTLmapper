# input method

import os
import sys
import pdb

class input_obj:

    def __init__(self,filename):
        self.species_tree_filename = None
        self.ultra_tree_bool = False
        self.boots_file = None
        self.output_dir = None
        self.penalties_filename = None
        self.special = None
        self.outgroup = None
        self.luca = None
        self.hr_scaling = 0.0
        self.iterate = None
        self.penalty_weight = 0.0
        self.event_guide_fn = None

        # read input data
        self.ReadInput(filename)

    def ReadInput(self,filename):
        '''store arguments from input datafile'''
        input_file = open(filename,'r')
        for line in input_file:

            # skip comments
            if line[0] == "#":
                continue

            line = line.strip()
            arg_parts = line.split('=')
            code = arg_parts[0]
            arg = arg_parts[1]
            if code == 'species':
                self.species_tree_filename = arg
            elif code == 'gene':
                self.boots_file = arg
            elif code == 'output':
                self.output_dir = arg
                if os.path.isdir(self.output_dir):
                    print "Output directory exists.  Overwrite?"
                    overwrite_s = sys.stdin.readline().strip().lower()
                    if overwrite_s == 'y' or overwrite_s == 'yes':
                        print "Overwriting"
                        continue
                    else:
                        print "Not overwriting.  Quitting."
                        sys.exit(1)
                    #endif
                #endif
            elif code == 'event_guide':
                if arg == 'None':
                    self.event_guide_fn = None
                else:
                    self.event_guide_fn = arg
            elif code == 'penalties':
                self.penalties_filename = arg
            elif code == 'ultrametric':
                self.ultra_tree_bool = eval(arg)
            elif code == 'special':
                self.special = arg
            elif code == 'outgroup':
                if arg == "None":
                    self.outgroup = None
                else:
                    self.outgroup = arg
            elif code == 'luca':
                self.luca = arg
            elif code == 'hr_scaling':
                self.hr_scaling = float(arg)
            elif code == 'iterate':
                if arg == "None":
                    self.iterate = None
                else:
                    self.iterate = float(arg)

        # make sure that if you're doing HR, you've got your outgroup
        if self.special == 'hr':
            if self.outgroup is None:
                print "if you're doing HR, you need an outgroup"
                sys.exit(1)

        input_file.close()

