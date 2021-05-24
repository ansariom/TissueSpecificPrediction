#!/usr/bin/env bash

tr '\n' '@' | tr '>' '\n' | sort | uniq | tr '\n' '>' | tr '@' '\n' | grep -v -E '^>$'
