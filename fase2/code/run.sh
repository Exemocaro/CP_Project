#!/bin/bash

sudo make clean
sudo make all
sudo make runpar
cd output
diff cp_average.txt default_average.txt
diff cp_output.txt default_output.txt
cd ..