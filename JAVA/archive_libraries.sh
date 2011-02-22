#/bin/sh
CUR=$PWD

cd access
ant jar
mv bbcf_access.jar ../jar
cd ..

cd parser
ant jar
mv bbcf_parser.jar ../jar
cd ..

cd sqlite
ant jar
mv bbcf_sqlite.jar ../jar
cd ..

