######################## prepare ##########################
mkdir ml2d_result
mkdir ml2d_result/fc
mkdir ml2d_result/fc/data3
mkdir ml2d_result/fc/data4
mkdir ml2d_result/fc/data5
mkdir ml2d_result/fc/data6
#running ml2d
input_fn="../dataset/GTM_paper/fc/data3"
output_fn="./ml2d_result/fc/data3/data3_it050_ml2d_classes100"
nr_classes=100
pixel_size=2.0
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
nr_pool=50
mpirun -n 30 -f all_machines_ib -perhost 1  ./bin/ml2d -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter -pool $nr_pool > data3_ml2d.r

input_fn="../dataset/GTM_paper/fc/data4"
output_fn="./ml2d_result/fc/data4/data4_it050_ml2d_classes100"
nr_classes=100
pixel_size=2.0
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
nr_pool=50
mpirun -n 30 -f all_machines_ib -perhost 1  ./bin/ml2d -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter -pool $nr_pool > data4_ml2d.r

input_fn="../dataset/GTM_paper/fc/data5"
output_fn="./ml2d_result/fc/data5/data5_it050_ml2d_classes100"
nr_classes=100
pixel_size=2.0
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
nr_pool=50

mpirun -n 30 -f all_machines_ib -perhost 1  ./bin/ml2d -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter -pool $nr_pool > data5_ml2d.r
input_fn="../dataset/GTM_paper/fc/data6"
output_fn="./ml2d_result/fc/data6/data6_it050_ml2d_classes100"
nr_classes=100
pixel_size=2.0
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
nr_pool=50
mpirun -n 30 -f all_machines_ib -perhost 1  ./bin/ml2d -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter -pool $nr_pool > data6_ml2d.r

######################## prepare ##########################
cp ./ml2d_result/fc/data3/data3_it050_ml2d_classes100.star ../dataset/GTM_paper/fc/
cp ./ml2d_result/fc/data3/data4_it050_ml2d_classes100.star ../dataset/GTM_paper/fc/
cp ./ml2d_result/fc/data3/data5_it050_ml2d_classes100.star ../dataset/GTM_paper/fc/
cp ./ml2d_result/fc/data3/data6_it050_ml2d_classes100.star ../dataset/GTM_paper/fc/
mkdir gtm_result
mkdir gtm_result/fc
mkdir gtm_result/fc/data3
mkdir gtm_result/fc/data4
mkdir gtm_result/fc/data5
mkdir gtm_result/fc/data6
#running gtm
input_fn="../dataset/GTM_paper/fc/data3_it050_ml2d_classes100"
output_fn="./gtm_result/fc/data3/data3_it050_gtm_classes100"
nr_classes=100
pixel_size=2
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
mpirun -n 10 -f all_machines_phi -perhost 1  ./bin/gtm -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter > data3_gtm.r

input_fn="../dataset/GTM_paper/fc/data4_it050_ml2d_classes100"
output_fn="./gtm_result/fc/data4/data4_it050_gtm_classes100"
nr_classes=100
pixel_size=2
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
mpirun -n 10 -f all_machines_phi -perhost 1  ./bin/gtm -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter > data4_gtm.r

input_fn="../dataset/GTM_paper/fc/data5_it050_ml2d_classes100"
output_fn="./gtm_result/fc/data5/data5_it050_gtm_classes100"
nr_classes=100
pixel_size=2
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
mpirun -n 10 -f all_machines_phi -perhost 1  ./bin/gtm -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter > data5_gtm.r

input_fn="../dataset/GTM_paper/fc/data6_it050_ml2d_classes100"
output_fn="./gtm_result/fc/data6/data6_it050_gtm_classes100"
nr_classes=100
pixel_size=2
nr_iter=50
nr_images=`grep mrcs $input_fn'.star' | wc -l`
mpirun -n 10 -f all_machines_phi -perhost 1  ./bin/gtm -i $input_fn -o $output_fn -n $nr_images -k $nr_classes -pixel $pixel_size -iter $nr_iter > data6_gtm.r
