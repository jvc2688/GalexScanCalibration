name=AIS_GAL_SCAN_03245_0001
#dis=/scratch/dw1519/galex/data/distortion/5-185
dis=/scratch/dw1519/galex/data/distortion/194-437
#dis=/scratch/dw1519/galex/data/distortion/446-887
#dis=/scratch/dw1519/galex/data/distortion/896-1157
#dis=/scratch/dw1519/galex/data/distortion/1166-1373
out=/scratch/dw1519/galex/co/co3245_1-10
echo $name
$HOME/dw1519/anaconda/bin/python co_rel_rot_csv_new_sec_dis_val_scatter.py ${name}-cal-sec 1 0.5 0.5 ${dis}-xa_low-centroids.npy ${dis}-xa_high-centroids.npy 1 400 100 ${out}




