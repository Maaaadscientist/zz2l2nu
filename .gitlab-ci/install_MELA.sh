#!/bin/bash

git clone https://github.com/JHUGen/JHUGenMELA.git
cd JHUGenMELA
git checkout 00cc82efec77a8dbc7c908f4f8203e5693e20e97

./setup.sh -j $(nproc) standalone
