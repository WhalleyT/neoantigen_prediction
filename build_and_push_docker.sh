#!/usr/bin/env bash

docker build -t neoantigen .
docker tag neoantigen twhalley93/neoantigen
docker push twhalley93/neoantigen