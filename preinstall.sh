#!/bin/sh
mkdir bin
#download picard
wget https://github.com/broadinstitute/picard/releases/download/2.9.4/picard.jar
mv picard.jar ./bin
#download GATK
wget https://github.com/broadinstitute/gatk/releases/download/4.beta.1/gatk-4.beta.1.zip
unzip gatk-4.beta.1.zip
rm gatk-4.beta.1.zip
mv ./gatk-4.beta.1/gatk-package-4.beta.1-local.jar ./bin
rm -rf ./gatk-4.beta.1
#complile bwa
cd ./RequiredProgram/bwa
make
mv bwa ../../bin/
cd ../../
#download samtools
wget https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2
tar jxvf samtools-1.5.tar.bz2
cd ./samtools-1.5
./configure --without-curses --disable-lzma --prefix=`pwd`
make
make install
mv samtools ../bin/
cd ..
rm samtools-1.5.tar.bz2
rm -rf samtools-1.5
#download sbt
wget https://github.com/sbt/sbt/releases/download/v0.13.15/sbt-0.13.15.tgz
tar zxvf sbt-0.13.15.tgz
mv ./sbt/bin/* ./bin
rm sbt-0.13.15.tgz
rm -rf sbt
#compile fgbio
cd ./RequiredProgram/fgbio
../../bin/sbt assembly
cp ./target/scala-2.12/fgbio-0.2.0-SNAPSHOT.jar ../../bin
cd ../../
#compile HapCut2
cd ./RequiredProgram/HapCut2
make
cp ./build/extractHAIRS ../../bin/
cp ./build/HAPCUT2 ../../bin/
cd ../../
#dowload supernova
cd ./bin
wget -O supernova-1.2.0.tar.gz "http://cf.10xgenomics.com/releases/assembly/supernova-1.2.0.tar.gz?Expires=1499756314&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9hc3NlbWJseS9zdXBlcm5vdmEtMS4yLjAudGFyLmd6IiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNDk5NzU2MzE0fX19XX0_&Signature=Th04LKpgZMo1bD6yB9Hj0QnBxdG1CTBfewKSR02p5Ykj6MvyFAlh2UDRw1qyLBI5V4cK3N2DD~e3GgheJuyzExFcOaV9aB~pWZcu07t45gGvBzOrO8AGN3kwG43sA9cwIDzoTc40P1MaoD8xNgfk8CGtycbOSfPc3qn99ZY1P6QoN4STnBFsjbUdRZP~yWNMZunhiUPbFVD5OU4j~JVtGV~WQDpHkzrxbn8V4NrENzav-d1ZSaBWr3aslIBKtoHc1olB6CkztLG5DF5kcCygYt6HoVMrRMasg95S8LVHja-UR9mZz3BTOjHT1SpGfCYGHvRJhTOm2LiUUtXh3QsctQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
tar zxvf supernova-1.2.0.tar.gz 
rm supernova-1.2.0.tar.gz
cd ..
#SPAdes
cd ./bin
wget http://cab.spbu.ru/files/release3.10.1/SPAdes-3.10.1-Linux.tar.gz
tar zxvf SPAdes-3.10.1-Linux.tar.gz
cp ./SPAdes-3.10.1-Linux/bin/spades.py ./
rm SPAdes-3.10.1-Linux.tar.gz
rm -rf ./SPAdes-3.10.1-Linux/
cd ..
#Install boost
cd ./RequiredProgram
mkdir boostlib
cd boostlib
wget http://downloads.sourceforge.net/project/boost/boost/1.56.0/boost_1_56_0.tar.bz2
tar jxf boost_1_56_0.tar.bz2
rm boost_1_56_0.tar.bz2
cd ../../
#install sparsehash
cd ./RequiredProgram/sparsehash/
./configure --prefix=`pwd`
make
make install
cd ../../
#Abyss
cd ./RequiredProgram
PATHname1=`pwd`
cd ./abyss
sparsehash="-I"${PATHname1}"/sparsehash/include/"
./autogen.sh
./configure CPPFLAGS=$sparsehash --prefix=${PATHname1}"/abyss-bin" --with-boost=${PATHname1}"/boostlib/boost_1_56_0/boost" --disable-openmp
make 
make install
cp ../abyss-bin/bin/abyss-pe ../../bin
rm -rf ../abyss-bin
cd ../../
#ARCS
cd ./RequiredProgram
PATHname1=`pwd`
cd ./ARCS
./autogen.sh 
./configure --with-boost=${PATHname1}"/boostlib/boost_1_56_0" --prefix=${PATHname1}"/ARCS-bin"
make install
cp ../ARCS-bin/bin/arcs ../../bin 
rm -rf ../ARCS-bin
cd ../../
##LINKS
wget http://www.bcgsc.ca/platform/bioinfo/software/links/releases/1.8.5/links_v1-8-5.tar.gz
tar zxvf links_v1-8-5.tar.gz
rm links_v1-8-5.tar.gz
cp ./links_v1.8.5/LINKS.pl ./bin
rm -rf ./links_v1.8.5
##Quast
cd ./RequiredProgram/quast
cp quast.py ../../bin
