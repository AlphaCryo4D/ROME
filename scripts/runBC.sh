######################## prepare ##########################
mkdir ml2d_result
mkdir ml2d_result/BC
mkdir ml2d_result/BC/class1
mkdir ml2d_result/BC/class2
mkdir ml2d_result/BC/class3
mkdir ml2d_result/BC/class5
mkdir ml2d_result/BC/class6
#running ml2d
input_fn="../dataset/GTM_paper/BC/class1_ori"
output_fn="./ml2d_result/BC/class1/class1_it050_ml2d_classes100"
nr_classes=100
pixel_size=1.74
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
nr_pool=50
mpirun -n 30 -f all_machines_ib -perhost 1  ./bin/ml2d -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter -pool $nr_pool > class1_ml2d.r

input_fn="../dataset/GTM_paper/BC/class2_ori"
output_fn="./ml2d_result/BC/class2/class2_it050_ml2d_classes100"
nr_classes=100
pixel_size=1.74
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
nr_pool=50
mpirun -n 30 -f all_machines_ib -perhost 1  ./bin/ml2d -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter -pool $nr_pool > class2_ml2d.r

input_fn="../dataset/GTM_paper/BC/class3_ori"
output_fn="./ml2d_result/BC/class3/class3_it050_ml2d_classes100"
nr_classes=100
pixel_size=2.0
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
nr_pool=50
mpirun -n 30 -f all_machines_ib -perhost 1  ./bin/ml2d -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter -pool $nr_pool > class3_ml2d.r

input_fn="../dataset/GTM_paper/BC/class5_ori"
output_fn="./ml2d_result/BC/class5/class5_it050_ml2d_classes100"
nr_classes=100
pixel_size=1.74
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
nr_pool=50
mpirun -n 30 -f all_machines_ib -perhost 1  ./bin/ml2d -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter -pool $nr_pool > class5_ml2d.r

input_fn="../dataset/GTM_paper/BC/class6_ori"
output_fn="./ml2d_result/BC/class6/class6_it050_ml2d_classes100"
nr_classes=100
pixel_size=1.74
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
nr_pool=50
mpirun -n 30 -f all_machines_ib -perhost 1  ./bin/ml2d -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter -pool $nr_pool > class6_ml2d.r
######################## prepare ##########################
cp ./ml2d_result/BC/class1/class1_it050_ml2d_classes100.star ../dataset/GTM_paper/BC/
cp ./ml2d_result/BC/class2/class2_it050_ml2d_classes100.star ../dataset/GTM_paper/BC/
cp ./ml2d_result/BC/class3/class3_it050_ml2d_classes100.star ../dataset/GTM_paper/BC/
cp ./ml2d_result/BC/class5/class5_it050_ml2d_classes100.star ../dataset/GTM_paper/BC/
cp ./ml2d_result/BC/class6/class6_it050_ml2d_classes100.star ../dataset/GTM_paper/BC/
mkdir gtm_result
mkdir gtm_result/BC
mkdir gtm_result/BC/class1
mkdir gtm_result/BC/class2
mkdir gtm_result/BC/class3
mkdir gtm_result/BC/class5
mkdir gtm_result/BC/class6
#running gtm
input_fn="../dataset/GTM_paper/BC/class1_it050_ml2d_classes100"
output_fn="./gtm_result/BC/class1/class1_it050_gtm_classes100"
nr_classes=100
pixel_size=1.74
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
mpirun -n 10 -f all_machines_phi -perhost 1  ./bin/gtm -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter > class1_gtm.r

input_fn="../dataset/GTM_paper/BC/class2_it050_ml2d_classes100"
output_fn="./gtm_result/BC/class2/class2_it050_gtm_classes100"
nr_classes=100
pixel_size=1.74
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
mpirun -n 10 -f all_machines_phi -perhost 1  ./bin/gtm -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter > class2_gtm.r

input_fn="../dataset/GTM_paper/BC/class3_it050_ml2d_classes100"
output_fn="./gtm_result/BC/class3/class3_it050_gtm_classes100"
nr_classes=100
pixel_size=1.74
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
mpirun -n 10 -f all_machines_phi -perhost 1  ./bin/gtm -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter > class3_gtm.r

input_fn="../dataset/GTM_paper/BC/class5_it050_ml2d_classes100"
output_fn="./gtm_result/BC/class5/class5_it050_gtm_classes100"
nr_classes=100
pixel_size=1.74
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
mpirun -n 10 -f all_machines_phi -perhost 1  ./bin/gtm -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter > class5_gtm.r

input_fn="../dataset/GTM_paper/BC/class6_it050_ml2d_classes100"
output_fn="./gtm_result/BC/class6/class6_it050_gtm_classes100"
nr_classes=100
pixel_size=1.74
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
mpirun -n 10 -f all_machines_phi -perhost 1  ./bin/gtm -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter > class6_gtm.r