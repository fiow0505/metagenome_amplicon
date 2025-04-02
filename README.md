# metagenome_amplicon (V3V4)

## Setting the analysis server (Linux server)
https://www.anaconda.com/docs/getting-started/miniconda/install
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh
conda env create –file metagenome_amplicon/metagenome_amplicon.yaml

conda activate metagenome_amplicon
cp -r metagenome_amplicon/script/ ${CONDA_PREFIX}/opt/metagenome_amplicon/
chmod +x ${CONDA_PREFIX}/opt/metagenome_amplicon/*
cp -r ${CONDA_PREFIX}/opt/metagenome_amplicon/*.sh ${CONDA_PREFIX}/bin/

git clone https://github.com/Zymo-Research/figaro.git
mv figaro ${CONDA_PREFIX}/opt/
pip3 install -r ${CONDA_PREFIX}/opt/figaro/requirements.txt
```

## install R packages
```
R
install.packages(c("cowplot", "ggprism", "ggpubr"))
```
