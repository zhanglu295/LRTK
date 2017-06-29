#!/bin/sh
mkdir bin
#download picard
wget https://github.com/broadinstitute/picard/releases/download/2.9.4/picard.jar
mv picard.jar ./bin
#download sbt
wget https://github.com/sbt/sbt/releases/download/v0.13.15/sbt-0.13.15.tgz
tar zxvf sbt-0.13.15.tgz
mv ./sbt/bin/* ./bin
rm sbt-0.13.15.tgz
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
