#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import random
import copy
from InitTest import *
from Truss.ClampedBeam import Categ, nEls

def neighbors(bb_args):
    float_args = []
    output = []
    for k in range(nEls):
        C = copy.copy(Categ)
        float_args.append(float(bb_args[k]))
        random.choice(C.remove(float_args[k]))
        output.append(random.choice(Categ))
    return (output)

def obj_func(bb_args):
	float_args = []
	for k in range(len(bb_args)):
		float_args.append(float(bb_args[k]))
	output = objectives_categ(float_args)
	return (output)

if __name__ == "__main__":
	input_args_txt = sys.argv[1]
	args = []
	for line in open(input_args_txt,"r"):
		for v in line.split(" "):
			args.append(v)
	outputs_bb = obj_func(args)
	buffer = "\n".join([str(o) for o in outputs_bb]+[" "])
	sys.stdout.write(buffer)
	sys.exit(0)
    