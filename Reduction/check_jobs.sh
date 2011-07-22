#!/bin/bash

tail -n 100 $1/*.log | egrep -B 2 '=>|rfcp|Red|Break' | egrep  '=>|rfcp|[Rr]ed|Break'

