#!/bin/sh
home=`pwd -P`;
#Change the following line to the path to the new pool!
for i in `find * -type l`; 
do target=$(readlink -- "${i}"); 
target=${target/*\/StatsProject\/16S/\/gpfs1\/projects\/Pools\/16SV4};
unlink ${i};
ln -s "${target}" "${i}"; 
done;
cd ${home};
exit 0;
