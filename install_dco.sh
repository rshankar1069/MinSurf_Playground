#! /bin/bash

echo 'DCL6I34NG ANNUAL CR11550:10927342 2020/03/31 "vwjQBc9GnzNfvHnfpc+bjojeaZ"' > nag_key.txt

# download and install dco
wget https://nag.co.uk/downloads/impl/dcl6i34ngl.tgz
tar -xzf dcl6i34ngl.tgz
./dcl6i34ngl/install.sh -silent -accept -installdir=dco_cpp

mv dco_cpp/dcl6i34ngl/* dco_cpp
rm -rf dco_cpp/dcl6i34ngl

# remove installer
rm -rf dcl6i34ngl dcl6i34ngl.tgz

