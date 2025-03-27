# install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
# remove sh file 
rm Miniconda3-latest-Linux-x86_64.sh

# create the environment 
conda create --name otters python=3.9 pandas=1.4.4 numpy=1.21.5 scipy=1.7.3 pip
# deactivate the conda environment
conda deactivate

 # activate the environment
conda activate otters

# install pysam
pip install pysam

# deactivate the environment
conda deactivate