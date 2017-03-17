if [ -z "$1" ]; then
    echo "Please supply original target"
    exit -1
fi

if [ -z "$2" ]; then
    echo "Please supply patch dir"
    exit -1
fi

if [ ! -d "$1" ]; then
    echo "$1 not found"
    exit -1
fi

if [ ! -d "$2" ]; then
    echo "$2 not found"
    exit -1
fi

if [ -d src ]; then
    echo "src already exists"
    exit -1
fi

ln -s $1/doc
ln -s $1/bindings
if [ -e $1/fallbacks ]; then 
  ln -s $1/fallbacks
fi
if [ -e $1/plugins ]; then
  ln -s $1/plugins
fi
ln -s $1/util
ln -s $1/tests
cp -r $1/config .
chmod -R u+w  config
if [ -e $1/prereqs ]; then
    ln -s $1/prereqs
fi

cp $1/Makefile* .
cp $1/config?* .
cp $1/aclocal.m4 .

chmod u+w Makefile*
chmod u+w config?*

mkdir src
cd src

mkdir mods
mkdir incs
mkdir libs

cp $1/src/mods/* mods
cp $1/src/incs/* incs
cp $1/src/libs/* libs

cp $1/src/Makefile* .

for l in $1/src/??_*; do
    DIRNAME=$(echo $l | sed s:$1/src/:: )
    if [ ! -e $2/$DIRNAME ]; then
	ln -s $l
    else
	mkdir $DIRNAME
	cd $DIRNAME
	for k in $l/*.F90; do
	    cp $l/Makefile* .
	    cp $l/abini?.in . 2> /dev/null
	    FILENAME=$(echo $k | sed s:$l/:: )
	    if [ ! -e $2/$DIRNAME/$FILENAME ]; then
		ln -s $k
		cp $(echo $k | sed s/F90/o/ ) . 2> /dev/null	
	    fi
	done
	chmod u+w *.o
	if [ -n "$(ls $2/$DIRNAME)" ]; then
	    MFILES=""
	    if [ -n "$(ls $2/$DIRNAME/*.f 2> /dev/null)" ]; then
		MFILES+=$(ls $2/$DIRNAME/*.f)" "
	    fi
	    if [ -n "$(ls $2/$DIRNAME/*.F90 2> /dev/null)" ]; then
		MFILES+=$(ls $2/$DIRNAME/*.F90)" "
	    fi
	    for k in $MFILES; do
		FILENAME=$(echo $k | sed s:$2/$DIRNAME/:: )
#		cp $k .
		ln -s $k

		LIBSNAME=$(grep -A 1 "Regular source files" Makefile.am | tail -n 1 | awk '{print $1}')

		OBJNAME=$(echo $FILENAME | sed s:F90:'$(OBJEXT)': | sed s:f$:'$(OBJEXT)':)
		if [ -z "$(grep $FILENAME Makefile)" ]; then
		    cat Makefile.in | sed s/"^am__objects_1 ="/"am__objects_1 = $OBJNAME"/ > tmp
		    mv tmp Makefile.in

		    cat Makefile.am | sed s/"^$LIBSNAME ="/"$LIBSNAME = $FILENAME"/ > tmp
		    mv tmp Makefile.am
		fi
	    done
	    cat Makefile.in | sed s/"^.SUFFIXES: .F90"/".SUFFIXES: .f .F90"/ \
		            | sed s/"^.F90.o:"/".f.o:\\n	\$(FC) -c -o \$\@ \$<\\n.F90.o:"/\
                    > tmp
	    mv tmp Makefile.in
	    echo "checking for $2/$DIRNAME/depends"
	    if [ -e "$2/$DIRNAME/depends" ]; then
		echo "found checking for $2/$DIRNAME/depends"
		cat $2/$DIRNAME/depends >> Makefile.in
		echo " " >> Makefile.in
		cat $2/$DIRNAME/depends >> Makefile.am
		echo " " >> Makefile.am
	    fi
	fi
	cd ..
    fi	
done
    
cd ..
