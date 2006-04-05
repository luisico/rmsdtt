#!/bin/sh

version=$1
echo "Packing version ${version:?}"

plugin=rmsdtt
dir=$plugin
tar=$plugin-v$version.tgz

# re-arrange documentation
mkdir doc
cp index.html rmsdtt_2.0.png doc
wget -q -O doc/vmdbackup.css http://physiology.med.cornell.edu/resources/physbio.css http://physiology.med.cornell.edu/faculty/hweinstein/vmdplugins/vmd.css
chmod -R ugo+rX doc

# Compress
cd ../
tar zcvf $tar $dir/pkgIndex.tcl $dir/rmsdtt.tcl $dir/doc
mv $tar $dir/versions
chmod 644 $dir/versions/$tar
cd $dir

rm -rf doc

