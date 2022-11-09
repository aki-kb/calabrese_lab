# calabrese_lab
Scripts and commands used for Braceros et al. 2022 manuscript.

Hi-C analysis
Juicer
$ juicer.sh -y mm9_Arima.txt -C 90000000 -L 15840 -Q 15840

Mega
$ mega.sh -g mm9 -L 15840 -Q 15840

Viewpoint
java -Xmx4048m -jar juicer_tools.jar pre merged_nodups.txt output.hic mm9
