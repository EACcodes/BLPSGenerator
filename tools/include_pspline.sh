if [ -z "$1" ]; then
    echo "Please supply patch dir"
    exit -1
fi

if [ -z "$2" ]; then
    echo "Please supply pspline dir"
    exit -1
fi

if [ ! -e "$2/libpsplinetotal.a" ]; then
    echo "pspline library not found at $2"
    exit -1
fi

echo
echo Adding pspline from $2 
echo to abinit build $1
echo

cd $1

ln -s $2
cp $2/pspline/LINUX/mod/* src/mods

cat src/98_main/Makefile.in | sed 's:lib_netcdf_libs@:lib_netcdf_libs@ -L$(abinit_builddir)/pspline/ -lpsplinetotal:' > tmp
mv tmp src/98_main/Makefile.in

echo installed PSPLINE
