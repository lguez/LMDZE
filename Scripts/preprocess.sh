# This is a script in Bash.

trap 'exit 1' ERR

##set -x

cd ~/Documents/Informatique_fonctionnement/Programs/LMDZ4_program/Pre-processed

for director in bibio dyn3d filtrez phylmd
  do
  cd $director
  rm -f *
  cd ..
done

cd ~/Documents/Informatique_fonctionnement/Programs/LMDZ4_program/libf

for director in bibio dyn3d filtrez phylmd
  do
  cd $director
  for filename in *.F?(90)
  do
    echo $filename
    suffix90=${filename##*.F}
    base=${filename%.*}
    g95 -E -DCPP_PHYS -DCPP_IOIPSL -I../grid -I../dyn3d -I../phylmd \
	-I/home/guez/netcdf-3.6.1/include $filename \
	> ../../Pre-processed/$director/$base.f$suffix90
  done
  cd ..
done
