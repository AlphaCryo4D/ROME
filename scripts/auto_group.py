'''
@auther yongbeima
@date 2016.08.16
@description : use to add group information to star file
@usage : python auto_group.py *.star nr_particles nr_micrographs nr_groups(smaller than nr_micrographs)
'''
import sys
import random
import math

filename = sys.argv[1]
particles_number = int(sys.argv[2])
micrograph_number = int(sys.argv[3])
group_number = int(sys.argv[4])

starfile = open(filename, 'r')
group_col_index = 0
write_head_complete = False

particles_micrograph_index = [random.randint(1,micrograph_number) for i in range(particles_number)]
micrograph_number_per_group = float(micrograph_number)/float(group_number)
# particles_micrograph_index =  sorted(particles_micrograph_index)
index = 0

for line in starfile :
	col = line.split()
	if (len(col) < 2): # read header's header
		print line,
	elif ("_rln") in line : # read header
		group_col_index = int(col[1][1:])
		print line,
	else : # read data
		if not write_head_complete:
			print "_rlnGroupName #"+str(group_col_index+1)
			print "_rlnMicrographName #"+str(group_col_index+2)
			write_head_complete = True
		for v in col:
			print v,
		print int(math.ceil(float(particles_micrograph_index[index]) / float(micrograph_number_per_group)))," ",particles_micrograph_index[index]
		index += 1

starfile.close()