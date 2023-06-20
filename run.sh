#!/bin/bash
#g++ QRSDetecorOnline.cpp -o output && ./output
if [[ "$OSTYPE" == "msys"* ]]; then
    g++ QRSDetecorOnline.cpp -g -o build/output.exe && ./build/output.exe
else 
    g++ QRSDetecorOnline.cpp -o build/output && ./build/output
fi
 