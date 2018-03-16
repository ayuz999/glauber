#!/bin/bash

[ -d output ] && rm -f output/*
[ -d output ] || mkdir -p output

[ -d LOG ]  && rm -f LOG/*
[ -d LOG ] || mkdir -p LOG

[ -d LIST ]  && rm -f LIST/*
[ -d LIST ] || mkdir -p LIST

[ -d table ]  && rm -f table/*
[ -d table ] || mkdir -p table

[ -d figure ]  && rm -f figure/*
[ -d figure ] || mkdir -p figure
