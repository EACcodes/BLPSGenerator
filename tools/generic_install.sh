PATCH_DIR=$PWD

if [ -z "$1" ]; then
    echo "Please do not call this file directly: use the install.sh provided in the respective directories"
    exit -1
fi
if [ -z "$2" ]; then
    echo "Please do not call this file directly: use the install.sh provided in the respective directories"
    exit -1
fi

PROGRAM_NAME=$1
ABINIT_VERSION=$2

PATH_TO_FILE="path_to_$(echo $ABINIT_VERSION)"

if [ ! -e $PATH_TO_FILE ]; then
    echo "Path-file $PATH_TO_FILE not found"	
    echo "Please enter directory of compiled original $ABINIT_VERSION"
    read -e inputline
    ORIG_DIR="$inputline"

    if [ ! -e "$ORIG_DIR" ]; then
	echo $ORIG_DIR not found
	exit -1
    fi

    echo Thank you.
    echo $ORIG_DIR > $PATH_TO_FILE 
else
    ORIG_DIR=$(cat $PATH_TO_FILE)
fi

if [ ! -e path_to_install ]; then
    echo "Path-file path_to_install not found"	
    echo "Please enter target directory for install of $PROGRAM_NAME"
    echo "This directory should not exist yet!"	
    read -e inputline
    TRUNK_DIR="$inputline"
    if [ -e "$TRUNK_DIR" ]; then
	echo "ERROR: target directory $TRUNK_DIR already exists"
	exit -1
    fi
    echo Thank you.
    echo $TRUNK_DIR > path_to_install 
else
    TRUNK_DIR=$(cat path_to_install)
fi

if [ -n "$(grep NEEDS_PSPLINE $PATCH_DIR/install.sh)" ]; then
    if [ ! -e path_to_pspline ]; then
        echo "Path-file path_to_pspline not found"	
	echo "Please enter directory containing compiled pspline"
	read -e inputline
	PATH_PSPLINE="$inputline"
	if [ ! -e $PATH_PSPLINE/libpsplinetotal.a ]; then
	    echo "Error: pspline library not found at $PATH_PSPLINE"
	    exit -1
	fi
	echo $PATH_PSPLINE > path_to_pspline
    else
	PATH_PSPLINE=$(cat path_to_pspline)       
    fi
fi

echo
echo
echo Installing $PROGRAM_NAME
echo original source dir $ORIG_DIR
echo patch dir $PATCH_DIR
echo into $TRUNK_DIR
if [ -n "$(grep NEEDS_PSPLINE $PATCH_DIR/install.sh)" ]; then
    echo Including PSPLINE from $PATH_PSPLINE
fi
echo 
echo "proceed (y/n)? "
echo '(if any of the above is incorrect, answer no, delete the appropriate path_to file, and run install again)'
read inputline

what="$inputline"

if [ ! "${what}" = "y" ]; then
    echo Aborting
    exit -1
fi

if [ ! -e "$ORIG_DIR" ]; then
    echo Error: abinit not found in $ORIG_DIR
    exit -1
fi
if [ -n "$(grep NEEDS_PSPLINE $PATCH_DIR/install.sh)" ]; then
    if [ ! -e $PATH_PSPLINE/libpsplinetotal.a ]; then
	echo "Error: pspline library not found at $PATH_PSPLINE"
	exit -1
    fi
fi
if [ -e "$TRUNK_DIR" ]; then
    echo "ERROR: target directory $TRUNK_DIR already exists"
    exit -1
fi

if mkdir -p $TRUNK_DIR; then
    echo succesfully created dir
else
    echo Error: could not create $TRUNK_DIR? 
    exit -1
fi

cd $TRUNK_DIR
TRUNK_DIR=$PWD # as the previous TRUNK_DIR might have been a relative path
bash $PATCH_DIR/../tools/applypatch.sh $ORIG_DIR $PATCH_DIR
if [ -n "$(grep NEEDS_PSPLINE $PATCH_DIR/install.sh)" ]; then
    echo "CALLING INCLUDE PSPLINE"
    bash $PATCH_DIR/../tools/include_pspline.sh $TRUNK_DIR $PATH_PSPLINE 
fi

if [ -e $PATCH_DIR/configure.call ]; then
    echo "configure.call exists, and will be used. Please remove or edit to change"
else
    echo "Trying to get abinit configuration from config.log"
    echo $(cat config.log | grep ./configure | head -n 1 | sed 's:\$::') > $PATCH_DIR/configure.call
fi

echo " "
echo " "
echo " "
echo "Calling configure as:"
cat $PATCH_DIR/configure.call
echo "Perhaps there are quotes missing in this because abinit does not save these"
echo "If this is the case, please edit the configure.call file"
echo "and include the quotes there"

echo "proceed (y/n)? "
echo '(if the above is incorrect, answer no, edit the configure.call file, and run install again)'
read inputline

what="$inputline"

if [ ! "${what}" = "y" ]; then
    echo Aborting
    exit -1
fi
bash $PATCH_DIR/configure.call
make
cp $ORIG_DIR/configure .
if [ -n "$(grep NEEDS_PSPLINE $PATCH_DIR/install.sh)" ]; then
    echo "CALLING INCLUDE PSPLINE"
    bash $PATCH_DIR/../tools/include_pspline.sh $TRUNK_DIR $PATH_PSPLINE 
fi
# prevent make from trying to remake Makefiles
touch Makefile.in 
make
cd $PATCH_DIR
ln -s $TRUNK_DIR/src/98_main/abinit abinit_$(echo $PROGRAM_NAME)
