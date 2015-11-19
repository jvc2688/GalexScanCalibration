import centroid_plot
import aplpy

if __name__ == '__main__':
	i = 560
	filename = '../fits/co/co_map%d_%d_zoom_nocal.fits'%(i, i+1)
	centroid = centroid_plot.find_centroid(filename)
	print centroid
	fig = aplpy.FITSFigure(filename)
	fig.add_label(centroid[0], centroid[1], 'X', color='red')
	fig.show_grayscale(invert=True)
	fig.tick_labels.set_xformat('d.ddd')
	fig.tick_labels.set_yformat('d.ddd')
	fig.recenter(0., 0., radius=0.01)
	fig.add_grid()
	fig.save('../plots/co_map%d_%d.png'%(i, i+1))