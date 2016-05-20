#!/bin/bash

for i in {1..429};
do
./BSM_Analysis /eos/uscms/store/user/lchaparr/DoubleMuon/DoubleMuon_Run2015D_PromptReco-v4/160210_020957/0000/OutTree_$i.root Syncro.root
done
