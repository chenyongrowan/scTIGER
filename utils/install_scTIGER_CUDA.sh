#!/bin/bash

conda create -n scTIGER
conda activate scTIGER

conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
pip install numpy pandas==1.4.3 matplotlib networkx argparse scipy scanpy argparse leidenalg
bambi==0.9.0 arviz numba==0.56.4

chmod +x run_scTIGER.py
unzip ./Data/ProstateCancer/Patient4_Benign_endothelial.zip -d ./Data/ProstateCancer
