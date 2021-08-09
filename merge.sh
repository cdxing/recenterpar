#!/bin/bash
date

cd production/
cp ../myhadd.cpp .
cp ../myhadd.C .

rm mer.csh
touch mer.csh
chmod +x mer.csh

path=$PWD
echo "The path is $path"

echo -n "root -b -q -x 'myhadd.cpp(" >> mer.csh
echo -n '"'$path'"' >> mer.csh
echo -n ")'" >> mer.csh

./mer.csh


