#!/bin/bash

libraries="$(basename libraries.*.*.tar.gz)"
tar xzf $libraries

cd libraries && sed -i '/# X11 library location/,/# GAMIT size dependent variables/c\
# Specific to Linux Ubuntu and Debian\
X11LIBPATH /usr/lib/x86_64-linux-gnu\
X11INCPATH /usr/include\
# GAMIT size dependent variables' Makefile.config

cd .. && tar czf $libraries libraries
rm -rf libraries
sed -i "s/set ans = $</set ans = 'y'/g" install_software
sed -i '$d' install_software
chmod +x install_software install_updates