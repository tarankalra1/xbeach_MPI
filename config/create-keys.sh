#!/bin/bash
echo "start"
ssh-keygen -t rsa
echo "key created"
ssh-copy-id -i ~/.ssh/id_rsa.pub xbeach@$CLUSTER_HOSTNAME
echo "copy finished"