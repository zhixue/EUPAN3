TOOL_ROOT=`dirname $0`

# quast
git clone https://github.com/ablab/quast
cd quast && ./setup.py install && cd $TOOL_ROOT

# minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make && cd $TOOL_ROOT

# cdhit
git clone https://github.com/weizhongli/cdhit
cd cdhit && make && cd $TOOL_ROOT

# gclust
git clone https://github.com/niu-lab/gclust
cd gclust && make && cd $TOOL_ROOT

# blast LATEST=2.11.0
wget -c https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz
tar xzvf ncbi-blast-2.11.0+-x64-linux.tar.gz






