#!/bin/bash
# Compile FloeDyn
# cd floedyn
python3 ./waf configure --gcc -o build_docker
python3 ./waf --target FLOE
# TODO compile other targets
# Run SSHD to allow connections
/usr/sbin/sshd -D

exec "$@"