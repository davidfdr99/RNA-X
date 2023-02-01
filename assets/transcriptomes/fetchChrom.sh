#!/bin/bash

# Script to fetch chromosome length for given genome

# UCSC Kent utilities fetchChromSizes required

fetchChromSizes ${1} | grep -v '_*_' > "${1}".chromSizes
