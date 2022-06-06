######################## prepare ##########################
mkdir ml2d_result
mkdir ml2d_result/Inf
#running ml2d
input_fn="../dataset/GTM_paper/Inf/ring11_all"
output_fn="./ml2d_result/Inf/ring11_it050_ml2d_classes50"
nr_classes=50
pixel_size=1.74
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
nr_pool=30
mpirun -n 10 -f all_machines_phi -perhost 1  ./bin/ml2d -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter -pool $nr_pool > ring11_all.r

######################## prepare ##########################
mkdir gtm_result
mkdir gtm_result/Inf
mkdir gtm_result/Inf/class8
mkdir gtm_result/Inf/class45
#view (drawimage.py) the result get by ML2D,then select the class8 and class45
awk '{if(NR <= 22 || $15 == 8) print $0}' ./ml2d_result/Inf/ring11_it050_ml2d_classes50.star > ./ml2d_result/Inf/ring11_class8.star
awk '{if(NR <= 22 || $15 == 45) print $0}' ./ml2d_result/Inf/ring11_it050_ml2d_classes50.star > ./ml2d_result/Inf/ring11_class45.star
#cp the star to mrcs data floder
cp ./ml2d_result/Inf/ring11_class8.star ../dataset/GTM_paper/Inf/
cp ./ml2d_result/Inf/ring11_class45.star ../dataset/GTM_paper/Inf/
#running the gtm for class8 and class45
input_fn="../dataset/GTM_paper/Inf/ring11_class8"
output_fn="./gtm_result/Inf/class8/ring11_class8_it050_gtm_classes100"
nr_classes=100
pixel_size=1.74
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
mpirun -n 3 -f all_machines_phi -perhost 1  ./bin/gtm -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter > ring11_class8.r

input_fn="../dataset/GTM_paper/Inf/ring11_class45"
output_fn="./gtm_result/Inf/class45/ring11_class45_it050_gtm_classes100"
nr_classes=100
pixel_size=1.74
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
mpirun -n 3 -f all_machines_phi -perhost 1  ./bin/gtm -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter > ring11_class45.r