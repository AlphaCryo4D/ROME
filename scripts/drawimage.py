'''
@auther yongbeima
@date 2015.09.13
@description : use to draw *.pmat,*.mrcs,*.dat file
@usage : python drawimage.py ./filename.mrcs
'''
import numpy
import matplotlib.pyplot as plot
import sys,os
import math
import struct
import time
from datetime import datetime

filename = sys.argv[1]

#this for dat
if(len(sys.argv) <= 2):
	ordering = 'ieee-be'
else:
	ordering = sys.argv[2]

(filename_base,filename_extension) = os.path.splitext(filename)
filename_base = filename_base.split('/')[-1]
path = os.path.dirname(filename)
print(path," ",filename_base," ",filename_extension)
timestr = "_"+datetime.now().strftime('%Y%m%d%H%M%S')

if(filename_extension == '.pmat'):
	pmatfile = open(filename, 'r')

	#read minimum net from file
	data = []
	for line in pmatfile:
	  col = line.split()
	  data += [float(i) for i in col]
	pmatfile.close()

	# print(data)
	# x = (np.random.rand(1754**2) < 0.5).astype(int)
	# im = data.reshape(N,N)
	# plot.gray()
	# plot.imshow(im);
	# #plot.imshow(x, cmap='gray', interpolation='nearest', vmin=0, vmax=255)
	# plot.savefig('/Users/yongbeima/Pictures/drawimages/textimage.png')
	# plot.show()

	# plot pmat data #
	N = math.sqrt(len(data))
	print(N)
	x = numpy.array(data)
	im = x.reshape(N,N)
	plot.gray()
	plot.imshow(im);
	#plot.imshow(x, cmap='gray', interpolation='nearest', vmin=0, vmax=255)
	# plot.savefig('/Users/yongbeima/Pictures/drawimages/'+filename_base+'.png')
	plot.savefig(path+'/'+filename_base+'.png')
	plot.show()

elif(filename_extension == ".mrcs"):
	mrcsfile = open(filename, 'rb')
	mrcsHead = mrcsfile.read(256*4)
	(nx,ny,N,mode) = struct.unpack("iiii", mrcsHead[:16])
	print(nx)
	print(ny)
	print(N)
	print(mode)
	# os.mkdir(r'/Users/yongbeima/Pictures/drawimages/'+filename_base)
	os.mkdir(r''+path+'/'+filename_base+timestr)
	for image_index in range(0,N):
		content = mrcsfile.read(4*nx*ny)
		data = struct.unpack("f"*nx*ny, content)
		x = numpy.array(data)
		im = x.reshape(nx,ny)
		plot.gray()
		plot.imshow(im);
		#plot.imshow(x, cmap='gray', interpolation='nearest', vmin=0, vmax=255)
		# plot.savefig('/Users/yongbeima/Pictures/drawimages/'+filename_base+'/'+filename_base+'_'+str(image_index)+'.png')
		plot.savefig(path+'/'+filename_base+timestr+'/'+filename_base+'_'+str(image_index+1)+'.png')
		print(sum(numpy.abs(x)))
		# plot.show()
	mrcsfile.close()

elif(filename_extension == ".dat"):
	mrcsfile = open(filename, 'rb')
	mrcsHead = mrcsfile.read(256*4)
	#get nz ny,assume ny == nx
	if(ordering == 'ieee-be'):
		datHead = struct.unpack(">"+"f"*256, mrcsHead)
	elif(ordering == 'ieee-le'):
		datHead = struct.unpack(""+"f"*256, mrcsHead)
	else:
		print("order is ieee-be(default) or ieee-le for dat file.")
		sys.exit(1)
	#the datHead[0] always equal to 1,and datHead[25] is number of images when using stack
	N = int(max(datHead[0],datHead[25]))
	ny = int(datHead[1])
	nx = int(datHead[11])
	mode = int(datHead[4])
	print(nx)
	print(ny)
	print(N)
	print(mode)
	#skip padding between head and data
	mrcsfile.read(int(datHead[21])-256*4)

	os.mkdir(r''+path+'/'+filename_base+timestr)
	for image_index in range(0,N):
		#when only single image stack!!!!!
		if(int(datHead[25]) >= 1):
			mrcsfile.read(int(datHead[21]))
		content = mrcsfile.read(4*nx*ny)
		if(ordering == 'ieee-be'):
			data = struct.unpack(">"+"f"*nx*ny, content)
		elif(ordering == 'ieee-le'):
			data = struct.unpack(""+"f"*nx*ny, content)
		else:
			print("order is ieee-be(default) or ieee-le for dat file.")
			sys.exit(1)

		x = numpy.array(data)
		im = x.reshape(nx,ny)
		plot.gray()
		plot.imshow(im);
		#plot.imshow(x, cmap='gray', interpolation='nearest', vmin=0, vmax=255)
		# plot.savefig('/Users/yongbeima/Pictures/drawimages/'+filename_base+'/'+filename_base+'_'+str(image_index)+'.png')
		plot.savefig(path+'/'+filename_base+timestr+'/'+filename_base+'_'+str(image_index+1)+'.png')
		print(sum(numpy.abs(x)))
		# plot.show()
	mrcsfile.close()

else:
	print('do nothing')












