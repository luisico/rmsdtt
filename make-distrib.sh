#!/bin/sh

version=$1
echo "Packing version ${version:?}"

plugin=rmsdtt
dir=$plugin
tar=$plugin-v$version.tgz

files=(
rmsdtt.tcl
#iterativefit.tcl
pkgIndex.tcl
doc
)

# re-arrange documentation
mkdir doc
sed -e "/html5-reset\.css/d" -e "/vmd\.css/d" index.html > doc/index.html
cp rmsdtt_3.0.png doc
cat ../html5-reset.css ../vmd.css > doc/$plugin.css
chmod -R ugo+rX doc

# Build list of files
for f in ${files[@]}; do
    if [ ! -r $f ]; then
	echo "ERROR: File \"$dir/$f\" is not readable!"
	exit 11
    fi
    dirfiles="$dirfiles $dir/$f"
done

# Compress
cd ../
tar zcvf $tar $dirfiles
mv $tar $dir/versions
chmod 644 $dir/versions/$tar
cd $dir

# Clean
rm -rf doc

