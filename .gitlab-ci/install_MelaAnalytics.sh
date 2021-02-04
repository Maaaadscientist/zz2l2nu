#!/bin/bash

git clone https://github.com/MELALabs/MelaAnalytics.git
cd MelaAnalytics
git checkout c81ac33828aa053228cc0ffa97a17ce6907402be

cd CandidateLOCaster
make -j $(nproc)
cd ../EventContainer
make -j $(nproc)
cd ../GenericMEComputer
make -j $(nproc)
