#!/bin/bash

input="input treemix file"
output="treemix output file"

treemix -i ${input} -root YRI -k 500 -bootstrap -o ${output}
