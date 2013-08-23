#!/bin/bash

# This script modifies BLAT headers to allow C++ linking to the BLAT library (jkOwnLib.a).

cd ../inc
for header_file in *.h
do
    echo $header_file
    mv $header_file $header_file.bak
    echo '#ifdef __cplusplus
extern "C" {
#endif' | cat - $header_file.bak > $header_file
    echo '#ifdef __cplusplus
}
#endif' >> $header_file
done
