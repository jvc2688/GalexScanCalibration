import pyfits
import math

def get_pos(file_name):
	hdu_list = pyfits.open(file_name)
	hdu0 = hdu_list[0]
	initial_pos, final_pos, center, sky_range = [], [], [], []
	initial_pos.append(hdu0.header['RA_1'])
	initial_pos.append(hdu0.header['DEC_1'])
	final_pos.append(hdu0.header['RA_20'])
	final_pos.append(hdu0.header['DEC_20'])
	return initial_pos, final_pos


def get_pos_range(file_name=None, name_list=None):
	names = []
	initial_pos, final_pos, center, sky_range = [], [], [], []
	if file_name == None and name_list != None:
		with open(name_list) as f:
			names = f.readlines()
	elif file_name != None:
		names.append(file_name)
	print names

	initial_pos, final_pos = get_pos(names[0])
	length = len(names)
	for i in range(1, length):
		file_name = names[i]
		tmp_initial, tmp_final = get_pos(file_name)
		for i in range(2):
			initial_pos[i] = tmp_initial[i] if math.fabs(final_pos[i] - tmp_initial[i])\
												> math.fabs(final_pos[i] - initial_pos[i])\
											else initial_pos[i]
			final_pos[i] = tmp_final[i] if math.fabs(tmp_final[i] - initial_pos[i])\
												> math.fabs(final_pos[i] - initial_pos[i])\
											else final_pos[i]
	
	center.append((initial_pos[0]+final_pos[0])/2.)
	center.append((initial_pos[1]+final_pos[1])/2.)
	sky_range.append(math.fabs(final_pos[0]-initial_pos[0]))
	sky_range.append(math.fabs(final_pos[1]-initial_pos[1]))
	print center, sky_range
	return center, sky_range
